"""Markdown report rendered from the evidence graph only. Every sentence and table row cites graph nodes
as [node-id]; `check.py` enforces it. Text coming from tools is sanitised so that it cannot open a new
sentence or a fake reference."""
from __future__ import annotations

from .rules import TIER_NAMES


def clean(s):
    return (str(s).replace("[", "(").replace("]", ")").replace("|", "/").replace(". ", "; ")
            .replace("\n", " ").strip().rstrip("."))


def num(x):
    return f"{x:,}"


def ev(x):
    return f"{x:.1e}" if x else "0"


def _seg(s):
    return f"{num(s[0])}-{num(s[1])}"


def sp_id(sseqid):
    """'sp|P12345|NAME_HUMAN' -> 'P12345 (NAME_HUMAN)'."""
    parts = sseqid.split("|")
    return clean(f"{parts[1]} ({parts[2]})" if len(parts) == 3 else sseqid)


def _sp_short(sp):
    h = sp["best_hit"]
    if not h:
        return "none"
    return f"{sp_id(h['sseqid'])} {h['pident']:.1f} % id / {h['qcovhsp']:.0f} % cov"


def candidate_section(i, n):
    """n: id -> node. Returns the lines of one candidate's section."""
    orf, tier, exp = n[f"orf:{i}"]["values"], n[f"tier:{i}"]["values"], n[f"exp:{i}"]["values"]
    sp, pf, rna = n[f"sprot:{i}"]["values"], n[f"pfam:{i}"]["values"], n[f"rna:{i}"]["values"]
    ctx = n.get(f"ctx:{i}", {}).get("values")
    para = n.get(f"para:{i}", {}).get("values")
    retro = n.get(f"retro:{i}", {}).get("values")
    segs = orf["cds_segments"]
    where = f"{orf['genbank_seq']}" + (f" (RefSeq {orf['refseq_seq']})" if orf["refseq_seq"] else "")
    out = [f"### {tier['rank']}. {i}: tier {tier['tier']} ({tier['name']})", "",
           f"Predicted locus on {where}, {_seg((orf['start'], orf['end']))} ({orf['strand']} strand), "
           f"{len(segs)} CDS segment{'s' if len(segs) != 1 else ''}, {orf['aa_len']} aa [orf:{i}].", ""]

    fors, against, weak_hit = [], [], None
    h = sp["best_hit"]
    if h:
        line = (f"Best Swiss-Prot hit {sp_id(h['sseqid'])}: {clean(h['protein_name'])} ({clean(h['organism']) or 'organism not given'}) at {h['pident']:.1f} % identity over "
                f"{h['qcovhsp']:.0f} % of the query, e-value {ev(h['evalue'])}, bitscore {h['bitscore']:.1f}")
        if sp["homolog"]:
            fors.append(line + f", which counts as a homolog under the fixed rule [sprot:{i}] [rule:homolog].")
        else:
            weak_hit = (line + f"; it is below the homolog rule's query coverage, so it is not counted as "
                        f"evidence [sprot:{i}] [rule:homolog].")
    if pf["domains"]:
        doms = ", ".join(f"{clean(d['acc'])} {clean(d['name'])} (aa {d['env_from']}-{d['env_to']}, "
                         f"score {d['score']:.1f})" for d in pf["domains"])
        fors.append(f"Pfam-A domains above the gathering threshold: {doms} [pfam:{i}].")
    if rna.get("supplied") and rna["overlapping_transcripts"]:
        fors.append(f"{len(rna['overlapping_transcripts'])} supplied transcript(s) overlap the locus on the "
                    f"same strand, matching {len(rna['junctions_matched'])} of {rna['predicted_junctions']} "
                    f"predicted junctions; RNA is recorded, not used for the tier [rna:{i}].")

    if ctx and ctx["class"] == "pseudogene":
        ps = [g for g in ctx["overlapping"] if "pseudo" in g["biotype"]]
        names = ", ".join(f"{clean(g['name']) or 'GeneID ' + g['gene_id']} ({clean(g['biotype'])}, "
                          f"{g['strand']} strand)" for g in ps)
        against.append(f"The reference annotates a pseudogene overlapping this locus: {names} [ctx:{i}].")
    if retro and retro["fired"]:
        p, par = para["best_hit"], para["parent"]
        against.append(f"Retrocopy signature: the candidate has a single CDS segment and its best paralog in "
                       f"the reference proteome, {clean(p['sseqid'])} (gene {clean(par['gene_name']) or par['gene_id']}, "
                       f"{par['cds_segments']} CDS segments, elsewhere in the genome), matches at "
                       f"{p['pident']:.1f} % identity over {p['qcovhsp']:.0f} % of the query "
                       f"[para:{i}] [retro:{i}] [rule:retrocopy].")
    paralog_note = None
    if para and para["best_hit"] and not (retro and retro["fired"]):
        p, par = para["best_hit"], para["parent"]
        gname = (clean(par["gene_name"]) or par["gene_id"]) if par else "unknown"
        paralog_note = (f"Closest paralog in the reference proteome is {clean(p['sseqid'])} (gene {gname}) at "
                        f"{p['pident']:.1f} % identity over {p['qcovhsp']:.0f} % of the query; the retrocopy "
                        f"rule did not fire, so it is not counted as evidence [para:{i}] [retro:{i}].")
    for fl in orf["flags"]:
        against.append(f"ORF flag: {fl} [orf:{i}].")

    out.append("**Evidence FOR**")
    out.append("")
    if fors:
        out.extend(f"- {x}" for x in fors)
    else:
        out.append(f"- No Swiss-Prot homolog and no Pfam domain were found [sprot:{i}] [pfam:{i}].")
    out.append("")
    out.append("**Evidence AGAINST**")
    out.append("")
    if against:
        out.extend(f"- {x}" for x in against)
    else:
        refs = " ".join(f"[{x}]" for x in (f"ctx:{i}", f"retro:{i}") if x in n) + f" [orf:{i}]"
        out.append(f"- No pseudogene, retrocopy signature or ORF flag was found {refs}.")
    out.append("")
    out.append("**Context**")
    out.append("")
    if ctx:
        others = [g for g in ctx["overlapping"] if "pseudo" not in g["biotype"]]
        if others:
            names = ", ".join(f"{clean(g['name']) or 'GeneID ' + g['gene_id']} ({clean(g['biotype'])}, "
                              f"{g['strand']} strand)" for g in others)
            out.append(f"- Other reference genes overlapping the locus span: {names} [ctx:{i}].")
        elif ctx["class"] == "nothing":
            out.append(f"- The reference annotates nothing at this locus, on either strand [ctx:{i}].")
        near = ctx["nearest_protein_coding"]
        if near:
            out.append(f"- Nearest reference protein-coding gene: {clean(near['name']) or 'GeneID ' + near['gene_id']} "
                       f"(GeneID {near['gene_id']}), {num(near['distance_bp'])} bp away [ctx:{i}].")
    else:
        out.append(f"- No reference annotation was supplied, so there is no locus context [orf:{i}] [rule:candidates].")
    if weak_hit:
        out.append(f"- {weak_hit}")
    if paralog_note:
        out.append(f"- {paralog_note}")
    out.append(f"- RNA evidence: {'supplied' if rna.get('supplied') else 'not supplied'} [rna:{i}].")
    ps = orf["predictor_own_score"]
    out.append(f"- The predictor's own score, not evidence: Carbon-A confidence {clean(ps['confidence']) or 'NA'}, "
               f"status {clean(ps['status']) or 'missing'} [orf:{i}].")
    out.append("")
    out.append("**Suggested validation**")
    out.append("")
    e = exp["experiment"]
    if exp["rule"] == 1:
        target = para["best_hit"]["sseqid"] if para and para["best_hit"] else "the closest paralog"
        txt = (f"{e[0].upper() + e[1:]}, with primers on positions where the candidate differs from "
               f"{clean(target)}")
        refs = f"[exp:{i}] [para:{i}]"
    elif exp["rule"] == 3:
        j = exp["junction"]
        txt = (f"{e}, with primers in the first two CDS segments in transcript order, {_seg(j['exon1'])} "
               f"and {_seg(j['exon2'])}; the first exon ends at {num(j['donor'])}, the second starts at "
               f"{num(j['acceptor'])}, and the predicted intron is {num(j['intron_bp'])} bp")
        txt = txt[0].upper() + txt[1:]
        refs = f"[exp:{i}] [orf:{i}]"
    else:
        txt = e[0].upper() + e[1:]
        refs = f"[exp:{i}] [orf:{i}]"
    out.append(f"- {txt}, chosen by experiment rule {exp['rule']} {refs} [rule:experiments].")
    reasons = "; ".join(clean(r) for r in tier["reasons"]) or "no evidence either way"
    out.append(f"- Tier {tier['tier']} ({tier['name']}) by the fixed tier rule: {reasons} [tier:{i}] [rule:tiers].")
    out.append("")
    return out


def render(graph, title):
    n = {x["id"]: x for x in graph["nodes"]}
    s = n["summary:run"]["values"]
    ref = "db:ref_gff" in n
    lines = [f"# gene-evidence report: {clean(title)}", ""]
    carbon = n["db:carbon_gff"]["database"][0]
    if ref:
        r = n["db:ref_gff"]["database"][0]
        lines.append(f"Carbon-A predictions for {clean(carbon['version'])} were compared with the reference "
                     f"annotation {clean(r['version'])} [db:carbon_gff] [db:ref_gff].")
        lines.append(f"Of {num(s['carbon_loci'])} Carbon-A loci, {num(s['dropped_no_identical_twin'])} sit on a "
                     f"sequence with no identical RefSeq twin, {num(s['overlap_reference_cds'])} overlap a reference "
                     f"protein-coding CDS on the same strand, and {num(s['candidates'])} are candidates "
                     f"[summary:run] [rule:candidates].")
    else:
        lines.append(f"Carbon-A predictions for {clean(carbon['version'])} were read with no reference annotation, "
                     f"so all {num(s['candidates'])} loci are candidates and there is no locus context "
                     f"[db:carbon_gff] [summary:run] [rule:candidates].")
    lines.append(f"This report shows the {num(s['shown'])} highest-ranked candidates [summary:run].")
    lines.append("Tiers come from a fixed rule and are a reading aid, not a claim: no ranking of such "
                 "candidates has been validated, and tier 1 is not claimed to be enriched for real genes "
                 "[rule:tiers].")
    tiers = ", ".join(f"tier {t} {TIER_NAMES[int(t)]} {num(c)}" for t, c in sorted(s["tiers"].items()))
    lines.append(f"Tier counts over all candidates: {tiers} [summary:run].")
    ctxs = ", ".join(f"{clean(k)} {num(v)}" for k, v in sorted(s["locus_context"].items()))
    lines.append(f"Locus context over all candidates: {ctxs}; retrocopy signatures "
                 f"{num(s['retrocopy_signatures'])} [summary:run].")
    lines += ["", "## Ranked summary", "",
              "| rank | locus | tier | locus context | Swiss-Prot best hit | Pfam families | experiment | nodes |",
              "|--:|---|---|---|---|---|---|---|"]
    for i in graph["ranking"]:
        t, sp, pf, ex = (n[f"{k}:{i}"]["values"] for k in ("tier", "sprot", "pfam", "exp"))
        ctx = n[f"ctx:{i}"]["values"]["class"] if f"ctx:{i}" in n else "no reference"
        if f"retro:{i}" in n and n[f"retro:{i}"]["values"]["fired"]:
            ctx += ", retrocopy"
        refs = " ".join(f"[{k}:{i}]" for k in ("tier", "ctx", "retro", "sprot", "pfam", "exp") if f"{k}:{i}" in n)
        lines.append(f"| {t['rank']} | {i} | {t['tier']} {t['name']} | {ctx} | {_sp_short(sp)} | "
                     f"{', '.join(pf['families']) or 'none'} | {ex['experiment']} | {refs} |")
    lines += ["", "## Candidates", ""]
    for i in graph["ranking"]:
        lines += candidate_section(i, n)
    return "\n".join(lines).rstrip("\n") + "\n"
