"""Unit tests for overlap, the retrocopy rule, tiering, experiments and the report-claim checker."""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import fixtures  # noqa: E402

from gene_evidence import check, io, rules, tools  # noqa: E402
from gene_evidence.candidates import select  # noqa: E402
from gene_evidence.cli import main  # noqa: E402
from gene_evidence.intervals import IntervalIndex, distance, union  # noqa: E402


# ---- overlap --------------------------------------------------------------------------------------

def test_union_merges_overlapping_and_touching():
    assert union([(10, 20), (1, 5), (6, 8), (15, 30)]) == [(1, 8), (10, 30)]


def test_overlap_bp():
    idx = IntervalIndex([(100, 200), (300, 400)])
    assert idx.overlap_bp([(50, 99)]) == 0
    assert idx.overlap_bp([(50, 100)]) == 1
    assert idx.overlap_bp([(150, 350)]) == 51 + 51
    assert idx.overlap_bp([(201, 299)]) == 0
    assert idx.overlap_bp([(150, 160), (155, 170)]) == 21  # query merged first


def test_distance():
    assert distance(1, 10, 5, 20) == 0
    assert distance(1, 10, 12, 20) == 1
    assert distance(30, 40, 1, 10) == 19


def test_candidate_selection(tmp_path):
    p = fixtures.write(str(tmp_path))
    loci = io.read_carbon_gff(p["carbon.gff3"])
    ref = io.read_reference_gff(p["ref.gff"])
    cands, counts = select(loci, ref, io.read_assembly_report(p["report.txt"]))
    ids = [c.locus.id for c in cands]
    # g7 overlaps a protein-coding CDS on its strand; g8 has no identical twin; g1's overlap is pseudo=true;
    # g9 overlaps a protein-coding CDS only on the other strand
    assert ids == ["g1", "g2", "g3", "g4", "g5", "g6", "g9"]
    assert counts == {"carbon_loci": 9, "dropped_no_identical_twin": 1, "overlap_reference_cds": 1,
                      "candidates": 7}
    assert {c.ref_seq for c in cands} == {"NC_1.1"}


def test_no_reference_means_every_locus_is_a_candidate(tmp_path):
    p = fixtures.write(str(tmp_path))
    cands, counts = select(io.read_carbon_gff(p["carbon.gff3"]))
    assert counts["candidates"] == 9 and counts["dropped_no_identical_twin"] == 0


def test_locus_context_classes(tmp_path):
    p = fixtures.write(str(tmp_path))
    genes = io.read_reference_gff(p["ref.gff"]).genes["NC_1.1"]
    assert rules.locus_context(1000, 1300, genes)["class"] == "pseudogene"
    assert rules.locus_context(50000, 50300, genes)["class"] == "non-coding gene"
    assert rules.locus_context(80000, 80400, genes)["class"] == "protein-coding gene"
    ctx = rules.locus_context(70000, 70150, genes)
    assert ctx["class"] == "nothing"
    assert ctx["nearest_protein_coding"] == {"gene_id": "15", "name": "Pc2", "strand": "+",
                                             "distance_bp": 79000 - 70150 - 1}


# ---- retrocopy rule ------------------------------------------------------------------------------

HIT = {"pident": 98.0, "qcovhsp": 95.0}


@pytest.mark.parametrize("segments,hit,parent_segs,overlaps,fired", [
    (1, HIT, 3, False, True),
    (2, HIT, 3, False, False),                                # candidate has an intron
    (1, {"pident": 89.9, "qcovhsp": 95.0}, 3, False, False),  # identity below 90 %
    (1, {"pident": 98.0, "qcovhsp": 79.0}, 3, False, False),  # coverage below 80 %
    (1, HIT, 1, False, False),                                # parent intronless too
    (1, HIT, 3, True, False),                                 # parent overlaps the candidate
    (1, None, None, None, False),                             # no paralog hit
])
def test_retrocopy_rule(segments, hit, parent_segs, overlaps, fired):
    assert rules.retrocopy(segments, hit, parent_segs, overlaps)["fired"] is fired


# ---- tiers and experiments ---------------------------------------------------------------------

@pytest.mark.parametrize("ctx,retro,homolog,fams,flags,tier", [
    ("pseudogene", False, True, 2, [], 5),     # AGAINST-strong overrides FOR
    ("nothing", True, True, 1, [], 5),
    ("nothing", False, True, 1, [], 1),
    ("nothing", False, True, 0, [], 2),
    ("non-coding gene", False, False, 1, ["no leading M"], 2),  # FOR beats a weak flag
    ("nothing", False, False, 0, ["internal stop"], 4),
    ("nothing", False, False, 0, [], 3),
    (None, False, False, 0, [], 3),            # no reference
])
def test_tiers(ctx, retro, homolog, fams, flags, tier):
    assert rules.tier(ctx, retro, homolog, fams, flags)[0] == tier


def test_homolog_needs_coverage():
    assert rules.is_homolog({"evalue": 1e-10, "qcovhsp": 50.0})
    assert not rules.is_homolog({"evalue": 1e-10, "qcovhsp": 49.9})
    assert not rules.is_homolog({"evalue": 1e-4, "qcovhsp": 90.0})


def test_orf_flags():
    assert rules.orf_flags("MAAA*", "complete") == []
    assert rules.orf_flags("MA*AA", "complete") == ["internal stop"]
    assert rules.orf_flags("AAA", "incomplete") == ["no leading M", "status incomplete"]


def test_experiment_rules():
    assert rules.experiment(5, True, HIT, [(1, 10)], "+", 50)["rule"] == 1
    assert rules.experiment(5, False, None, [(1, 10)], "+", 50)["rule"] == 2
    assert rules.experiment(2, False, None, [(1, 300)], "+", 99)["rule"] == 4
    assert rules.experiment(2, False, None, [(1, 300)], "+", 100)["rule"] == 5
    plus = rules.experiment(1, False, None, [(100, 200), (300, 400), (500, 600)], "+", 100)
    assert plus["rule"] == 3 and plus["junction"]["donor"] == 200 and plus["junction"]["acceptor"] == 300
    assert plus["junction"]["intron_bp"] == 99
    minus = rules.experiment(1, False, None, [(100, 200), (300, 400), (500, 600)], "-", 100)
    # transcript order on the minus strand starts at the highest coordinate
    assert minus["junction"]["donor"] == 500 and minus["junction"]["acceptor"] == 400
    assert minus["junction"]["intron_bp"] == 99


def test_parse_diamond_keeps_best_bitscore(tmp_path):
    f = tmp_path / "h.tsv"
    f.write_text(fixtures.SPROT_HITS)
    best = tools.parse_diamond(str(f))
    assert best["g3"]["sseqid"] == "sp|P3|D_HUMAN"
    assert tools.swissprot_title(best["g3"]["stitle"]) == ("Protein D [x]", "Homo sapiens", "9606")


# ---- report-claim checker ---------------------------------------------------------------------

NODES = ["orf:g1", "ctx:g1"]


def test_checker_passes_cited_text():
    text = "# Title\n\n**Evidence FOR**\n\n- One claim [orf:g1]. Another one [ctx:g1].\n\n" \
           "| a | b |\n|---|---|\n| x | [orf:g1] |\n"
    assert check.check(text, NODES) == []


def test_checker_fails_on_uncited_sentence():
    text = "- One claim [orf:g1]. This one has no reference.\n"
    problems = check.check(text, NODES)
    assert len(problems) == 1 and "uncited" in problems[0]


def test_checker_fails_on_dangling_reference():
    problems = check.check("A claim [orf:g2].\n", NODES)
    assert len(problems) == 1 and "not a graph node" in problems[0]


def test_checker_fails_on_uncited_table_row():
    problems = check.check("| a | b |\n|---|---|\n| x | y |\n", NODES)
    assert len(problems) == 1


# ---- end to end on synthetic fixtures (tools replaced by pre-seeded raw tables) -----------------

def _run(tmp_path, name, **kw):
    p = fixtures.write(str(tmp_path / "in"))
    out = str(tmp_path / name)
    fixtures.seed_raw(out, reference=kw.get("reference", True))
    rc = main(fixtures.argv(p, out, **kw))
    return rc, out


def test_end_to_end_tiers_and_traceability(tmp_path):
    rc, out = _run(tmp_path, "a")
    assert rc == 0
    g = json.load(open(os.path.join(out, "graph.json")))
    tiers = {n["values"]["locus"]: n["values"]["tier"] for n in g["nodes"] if n["type"] == "candidate"}
    assert tiers == {"g1": 5, "g2": 5, "g3": 1, "g4": 2, "g5": 4, "g6": 3, "g9": 4}
    assert g["ranking"] == ["g3", "g4", "g6", "g5", "g9", "g1", "g2"]
    for n in g["nodes"]:
        assert set(n) == {"id", "type", "tool", "database", "parameters", "values"}
    text = open(os.path.join(out, "report.md")).read()
    assert check.check(text, [n["id"] for n in g["nodes"]]) == []
    assert main(["check-report", os.path.join(out, "report.md"), os.path.join(out, "graph.json")]) == 0
    versions = json.load(open(os.path.join(out, "versions.json")))
    for k in ("carbon_gff", "carbon_proteins", "ref_gff", "ref_proteins", "swissprot", "pfam"):
        assert len(versions["inputs_and_databases"][k]["sha256"]) == 64


def test_end_to_end_is_byte_deterministic(tmp_path):
    _, a = _run(tmp_path, "a")
    _, b = _run(tmp_path, "b")
    for f in ("graph.json", "report.md"):
        assert open(os.path.join(a, f), "rb").read() == open(os.path.join(b, f), "rb").read()


def test_tampered_report_fails_the_checker(tmp_path):
    _, out = _run(tmp_path, "a")
    rp = os.path.join(out, "report.md")
    with open(rp, "a") as f:
        f.write("\nThis candidate is clearly a real gene.\n")
    assert main(["check-report", rp, os.path.join(out, "graph.json")]) == 1


def test_top_limits_report(tmp_path):
    _, out = _run(tmp_path, "a", top=2)
    g = json.load(open(os.path.join(out, "graph.json")))
    assert g["ranking"] == ["g3", "g4"]
    summary = next(n for n in g["nodes"] if n["id"] == "summary:run")["values"]
    assert summary["candidates"] == 7 and summary["shown"] == 2


def test_no_reference_run(tmp_path):
    rc, out = _run(tmp_path, "a", reference=False)
    assert rc == 0
    g = json.load(open(os.path.join(out, "graph.json")))
    assert not any(n["id"].startswith(("ctx:", "para:", "retro:")) for n in g["nodes"])
    assert "no reference annotation" in open(os.path.join(out, "report.md")).read()


def test_rna_absent_is_a_recorded_noop_and_present_matches_junctions(tmp_path):
    _, out = _run(tmp_path, "a")
    g = {n["id"]: n for n in json.load(open(os.path.join(out, "graph.json")))["nodes"]}
    assert g["rna:g3"]["values"] == {"supplied": False}
    _, out = _run(tmp_path, "b", rna=True)
    g = {n["id"]: n for n in json.load(open(os.path.join(out, "graph.json")))["nodes"]}
    assert g["rna:g3"]["values"]["junctions_matched"] == [[40501, 40799]]
    assert g["tier:g3"]["values"]["tier"] == 1  # RNA is recorded, not used for the tier
