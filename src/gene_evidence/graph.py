"""Evidence graph: one node per tool output or rule application, edges to what each node used.

Pure function of its inputs: no clock, no paths beyond file basenames, sorted keys, so two runs on the
same inputs serialise to identical bytes (acceptance criterion (a))."""
from __future__ import annotations

import json

from . import rna as rna_mod
from . import rules
from .tools import PARALOG_PARAMS, PFAM_PARAMS, SWISSPROT_PARAMS, swissprot_title

SCHEMA = "gene-evidence/evidence-graph/1"


def _node(nid, ntype, tool, database, parameters, values):
    return {"id": nid, "type": ntype, "tool": tool, "database": database,
            "parameters": parameters, "values": values}


def evaluate(cand, reference, protein, sprot_hit, paralog_hit, domains, rna_tx):
    """All evidence and rule outcomes for one candidate, as plain values."""
    L = cand.locus
    ctx = None
    if reference is not None:
        ctx = rules.locus_context(L.start, L.end, reference.genes.get(cand.ref_seq, []))
    parent = None
    if paralog_hit is not None and reference is not None:
        p = reference.proteins.get(paralog_hit["sseqid"])
        g = reference.gene_by_id.get(p["gene_id"]) if p else None
        parent = {"protein_id": paralog_hit["sseqid"], "gene_id": p["gene_id"] if p else "",
                  "gene_name": g.name if g else "", "cds_segments": p["segments"] if p else None,
                  "overlaps_candidate": (g.seq == cand.ref_seq and g.start <= L.end and L.start <= g.end)
                  if g else None}
    retro = None
    if reference is not None:
        retro = rules.retrocopy(len(L.cds), paralog_hit,
                                parent["cds_segments"] if parent else None,
                                parent["overlaps_candidate"] if parent else None)
    aa = protein.rstrip("*") if protein else ""
    flags = rules.orf_flags(protein, L.status)
    homolog = rules.is_homolog(sprot_hit)
    families = sorted({d["acc"] for d in domains or []})
    tier_no, reasons = rules.tier(ctx["class"] if ctx else None, bool(retro and retro["fired"]),
                                  homolog, len(families), flags)
    exp = rules.experiment(tier_no, bool(retro and retro["fired"]), paralog_hit, L.cds, L.strand, len(aa))
    rna = rna_mod.evidence(rna_tx, [L.seq, cand.ref_seq], L.strand, L.start, L.end, L.cds)
    return {"cand": cand, "ctx": ctx, "parent": parent, "retro": retro, "aa_len": len(aa),
            "starts_m": aa.startswith("M"), "internal_stop": "*" in aa, "flags": flags,
            "homolog": homolog, "families": families, "tier": tier_no, "reasons": reasons,
            "experiment": exp, "rna": rna, "sprot": sprot_hit, "paralog": paralog_hit,
            "domains": domains or [],
            "order": (tier_no, -len(families), -(sprot_hit["bitscore"] if sprot_hit else 0.0), L.id)}


def build(cands, counts, reference, proteins, sprot, paralog, pfam, rna_tx, dbs, tools, top=None):
    """dbs: key -> {"name", "file", "version", "sha256"} for every input file and database.
    tools: name -> version. sprot/paralog/pfam: parsed tables, or None if that stage did not run."""
    evals = [evaluate(c, reference, proteins.get(c.locus.id, ""),
                      sprot.get(c.locus.id) if sprot is not None else None,
                      paralog.get(c.locus.id) if paralog is not None else None,
                      pfam.get(c.locus.id, []) if pfam is not None else None, rna_tx)
             for c in cands]
    evals.sort(key=lambda e: e["order"])
    shown = evals if top is None else evals[:top]

    ge = {"name": "gene-evidence", "version": tools["gene-evidence"]}
    db = lambda *keys: [dbs[k] for k in keys if k in dbs]  # noqa: E731
    nodes, edges = [], []

    for k in sorted(dbs):
        nodes.append(_node(f"db:{k}", "input", None, [dbs[k]], {}, {}))
    nodes.append(_node("rule:candidates", "rule", ge, None, {}, {
        "text": "Carbon-A loci whose CDS overlaps no reference protein-coding CDS (pseudo=true excluded) "
                "on the same strand; GenBank->RefSeq names through assembly-report rows with "
                "relationship '='; without a reference every locus is a candidate",
        "reference_supplied": reference is not None}))
    nodes.append(_node("rule:retrocopy", "rule", ge, None, dict(rules.RETROCOPY), {}))
    nodes.append(_node("rule:homolog", "rule", ge, None, dict(rules.HOMOLOG), {}))
    nodes.append(_node("rule:tiers", "rule", ge, None, {}, {
        "rule": rules.TIER_RULE, "names": {str(k): v for k, v in rules.TIER_NAMES.items()},
        "status": "a fixed reading aid; no ranking of these candidates has been validated"}))
    nodes.append(_node("rule:experiments", "rule", ge, None, {"small_orf_aa": rules.SMALL_ORF_AA},
                       {"rule": rules.EXPERIMENT_RULE}))
    tier_counts = {str(t): sum(e["tier"] == t for e in evals) for t in rules.TIER_NAMES}
    ctx_counts = {}
    for e in evals:
        cls = e["ctx"]["class"] if e["ctx"] else "no reference"
        ctx_counts[cls] = ctx_counts.get(cls, 0) + 1
    nodes.append(_node("summary:run", "summary", ge, db(*sorted(dbs)), {"top": top}, {
        **counts, "tiers": tier_counts, "locus_context": ctx_counts,
        "retrocopy_signatures": sum(bool(e["retro"] and e["retro"]["fired"]) for e in evals),
        "shown": len(shown)}))

    for rank, e in enumerate(shown, 1):
        L = e["cand"].locus
        i = L.id
        ev = []

        def add(nid, ntype, tool, database, params, values, uses=()):
            nodes.append(_node(nid, ntype, tool, database, params, values))
            for u in uses:
                edges.append({"from": nid, "to": u, "rel": "uses"})
            ev.append(nid)

        add(f"orf:{i}", "orf", {"name": "Carbon-A (as published)", "version": dbs.get("carbon_gff", {}).get("version", "")},
            db("carbon_gff", "carbon_proteins"), {}, {
                "genbank_seq": L.seq, "refseq_seq": e["cand"].ref_seq, "strand": L.strand,
                "start": L.start, "end": L.end, "cds_segments": [list(s) for s in L.cds],
                "aa_len": e["aa_len"], "starts_m": e["starts_m"], "internal_stop": e["internal_stop"],
                "flags": e["flags"],
                "predictor_own_score": {"confidence": L.confidence, "status": L.status}})
        if e["ctx"] is not None:
            add(f"ctx:{i}", "locus_context", ge, db("ref_gff", "ref_report"),
                {"overlap": "candidate gene span, either strand",
                 "nearest": "reference protein-coding gene on the same sequence, either strand"}, e["ctx"])
            p = e["paralog"]
            add(f"para:{i}", "paralog_hit", {"name": "DIAMOND", "version": tools.get("diamond", "")},
                db("ref_proteins", "ref_gff"), dict(PARALOG_PARAMS),
                {"best_hit": {k: p[k] for k in ("sseqid", "pident", "qcovhsp", "evalue", "bitscore",
                                                  "length", "qlen", "slen")} if p else None,
                 "parent": e["parent"]})
            add(f"retro:{i}", "rule_application", ge, None, {}, e["retro"],
                uses=(f"orf:{i}", f"para:{i}", "rule:retrocopy"))
        h = e["sprot"]
        sp_vals = None
        if h:
            name, organism, taxid = swissprot_title(h["stitle"])
            sp_vals = {"sseqid": h["sseqid"], "protein_name": name, "organism": organism, "taxid": taxid,
                       **{k: h[k] for k in ("pident", "qcovhsp", "evalue", "bitscore", "length", "qlen",
                                            "slen")}}
        add(f"sprot:{i}", "homology", {"name": "DIAMOND", "version": tools.get("diamond", "")},
            db("swissprot"), dict(SWISSPROT_PARAMS),
            {"ran": sprot is not None, "best_hit": sp_vals, "homolog": e["homolog"]}, uses=("rule:homolog",))
        add(f"pfam:{i}", "domains", {"name": "pyhmmer", "version": tools.get("pyhmmer", "")},
            db("pfam"), dict(PFAM_PARAMS),
            {"ran": pfam is not None, "domains": e["domains"], "families": e["families"]})
        add(f"rna:{i}", "rna", ge, db("rna_gtf") if e["rna"]["supplied"] else None,
            {"input": "transcript GTF"}, e["rna"])
        add(f"tier:{i}", "rule_application", ge, None, {},
            {"tier": e["tier"], "name": rules.TIER_NAMES[e["tier"]], "reasons": e["reasons"], "rank": rank},
            uses=tuple(x for x in (f"ctx:{i}", f"retro:{i}", f"sprot:{i}", f"pfam:{i}", f"orf:{i}")
                       if x in ev) + ("rule:tiers",))
        add(f"exp:{i}", "rule_application", ge, None, {}, e["experiment"],
            uses=(f"tier:{i}", f"orf:{i}", "rule:experiments") + ((f"para:{i}",) if f"para:{i}" in ev else ()))
        nodes.append(_node(f"cand:{i}", "candidate", ge, None, {}, {"locus": i, "rank": rank, "tier": e["tier"]}))
        edges.extend({"from": f"cand:{i}", "to": n, "rel": "has_evidence"} for n in ev)
        edges.append({"from": f"cand:{i}", "to": "rule:candidates", "rel": "selected_by"})

    nodes.sort(key=lambda n: n["id"])
    edges.sort(key=lambda x: (x["from"], x["to"], x["rel"]))
    return {"schema": SCHEMA, "ranking": [e["cand"].locus.id for e in shown], "nodes": nodes, "edges": edges}


def dumps(graph):
    return json.dumps(graph, indent=1, sort_keys=True, ensure_ascii=True) + "\n"
