"""The fixed rules of stage 1 (BRIEF.md): locus context class, retrocopy signature, homolog, tiers,
validation experiment. Pure functions; every threshold lives in this file and is copied into the graph."""
from __future__ import annotations

from .intervals import distance, span_overlaps

RETROCOPY = {"candidate_max_cds_segments": 1, "min_identity_pct": 90.0, "min_query_coverage_pct": 80.0,
             "parent_min_cds_segments": 2, "parent_must_not_overlap_candidate": True}
HOMOLOG = {"max_evalue": 1e-5, "min_query_coverage_pct": 50.0}
SMALL_ORF_AA = 100

TIER_NAMES = {1: "FOR-strong", 2: "FOR", 3: "no evidence", 4: "AGAINST-weak", 5: "AGAINST-strong"}
TIER_RULE = [
    "5 AGAINST-strong: locus class is pseudogene, or the retrocopy signature fired (overrides FOR)",
    "1 FOR-strong: Swiss-Prot homolog and >= 1 Pfam domain",
    "2 FOR: Swiss-Prot homolog or >= 1 Pfam domain",
    "4 AGAINST-weak: any ORF flag and no FOR evidence",
    "3 no evidence: everything else",
    "order within a tier: distinct Pfam families desc, best Swiss-Prot bitscore desc, locus id asc",
]
EXPERIMENT_RULE = [
    "1 tier 5 with a retrocopy signature or a paralog hit: paralog-discriminating RT-PCR",
    "2 tier 5 without a paralog hit: long-read cDNA sequencing across the locus",
    "3 >= 2 CDS segments: RT-PCR across the first predicted splice junction",
    f"4 one CDS segment, protein < {SMALL_ORF_AA} aa: ribosome profiling",
    f"5 one CDS segment, protein >= {SMALL_ORF_AA} aa: long-read cDNA sequencing",
]


def locus_context(start, end, genes):
    """genes: RefGene list of the candidate's RefSeq sequence. Class precedence: pseudogene,
    protein-coding gene, non-coding gene, nothing."""
    overlapping = [g for g in genes if span_overlaps(start, end, g.start, g.end)]
    biotypes = {g.biotype for g in overlapping}
    if any("pseudo" in b for b in biotypes):
        cls = "pseudogene"
    elif "protein_coding" in biotypes:
        cls = "protein-coding gene"
    elif biotypes:
        cls = "non-coding gene"
    else:
        cls = "nothing"
    nearest = None
    for g in genes:
        if g.biotype != "protein_coding":
            continue
        d = distance(start, end, g.start, g.end)
        if nearest is None or (d, g.gene_id) < (nearest["distance_bp"], nearest["gene_id"]):
            nearest = {"gene_id": g.gene_id, "name": g.name, "strand": g.strand, "distance_bp": d}
    return {
        "class": cls,
        "overlapping": [{"gene_id": g.gene_id, "name": g.name, "biotype": g.biotype, "strand": g.strand,
                         "start": g.start, "end": g.end} for g in overlapping],
        "nearest_protein_coding": nearest,
    }


def retrocopy(n_segments, paralog, parent_segments, parent_overlaps):
    """paralog: best hit dict against the reference proteome, or None."""
    cond = {
        "single_cds_segment": n_segments <= RETROCOPY["candidate_max_cds_segments"],
        "paralog_hit": paralog is not None,
        "identity_ok": bool(paralog) and paralog["pident"] >= RETROCOPY["min_identity_pct"],
        "coverage_ok": bool(paralog) and paralog["qcovhsp"] >= RETROCOPY["min_query_coverage_pct"],
        "parent_multi_segment": parent_segments is not None
        and parent_segments >= RETROCOPY["parent_min_cds_segments"],
        "parent_elsewhere": parent_overlaps is False,
    }
    return {"fired": all(cond.values()), "conditions": cond}


def is_homolog(hit):
    return bool(hit) and hit["evalue"] <= HOMOLOG["max_evalue"] \
        and hit["qcovhsp"] >= HOMOLOG["min_query_coverage_pct"]


def orf_flags(protein, status):
    flags = []
    if not protein:
        flags.append("no protein sequence")
        return flags
    if not protein.startswith("M"):
        flags.append("no leading M")
    if "*" in protein.rstrip("*"):
        flags.append("internal stop")
    if status != "complete":
        flags.append(f"status {status or 'missing'}")
    return flags


def tier(context_class, retro_fired, homolog, n_families, flags):
    """-> (tier number, reasons). First match wins, in the order of TIER_RULE."""
    against = []
    if context_class == "pseudogene":
        against.append("annotated pseudogene at the locus")
    if retro_fired:
        against.append("retrocopy signature")
    if against:
        return 5, against
    fors = (["Swiss-Prot homolog"] if homolog else []) + (["Pfam domain"] if n_families else [])
    if len(fors) == 2:
        return 1, fors
    if fors:
        return 2, fors
    if flags:
        return 4, list(flags)
    return 3, []


def first_junction(cds, strand):
    """First splice junction in transcript order: donor/acceptor coordinates and intron length."""
    segs = sorted(cds)
    if strand == "-":
        a, b = segs[-1], segs[-2]
        donor, acceptor = a[0], b[1]
        return {"exon1": list(a), "exon2": list(b), "donor": donor, "acceptor": acceptor,
                "intron_bp": donor - acceptor - 1}
    a, b = segs[0], segs[1]
    return {"exon1": list(a), "exon2": list(b), "donor": a[1], "acceptor": b[0],
            "intron_bp": b[0] - a[1] - 1}


def experiment(tier_no, retro_fired, paralog, cds, strand, aa_len):
    if tier_no == 5 and (retro_fired or paralog):
        return {"rule": 1, "experiment": "paralog-discriminating RT-PCR"}
    if tier_no == 5:
        return {"rule": 2, "experiment": "long-read cDNA sequencing"}
    if len(cds) >= 2:
        return {"rule": 3, "experiment": "RT-PCR across the first predicted splice junction",
                "junction": first_junction(cds, strand)}
    if aa_len < SMALL_ORF_AA:
        return {"rule": 4, "experiment": "ribosome profiling"}
    return {"rule": 5, "experiment": "long-read cDNA sequencing"}
