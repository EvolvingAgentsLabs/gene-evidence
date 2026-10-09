"""`gene-evidence` command line: analyze, check-report, fetch."""
from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import sys
import time

from . import __version__, check, fetch, graph, io, report, rna, tools
from .candidates import select


def _rss_mb(who):
    r = resource.getrusage(who).ru_maxrss
    return round(r / (1024 * 1024) if sys.platform == "darwin" else r / 1024, 1)  # bytes on macOS, KiB on Linux


class Log:
    def __init__(self, path):
        self.f = open(path, "a", buffering=1)
        self.t0 = time.time()

    def __call__(self, msg):
        line = f"[{time.time() - self.t0:8.1f}s] {msg}"
        self.f.write(line + "\n")
        print(line, file=sys.stderr, flush=True)


def _db(name, path, version):
    return {"name": name, "file": os.path.basename(path), "version": version, "sha256": io.sha256(path)}


def analyze(a):
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    raw = os.path.join(a.out, "raw")
    os.makedirs(raw, exist_ok=True)
    log = Log(os.path.join(a.out, "gene-evidence.log"))
    stages = {}

    def stage(name, t):
        stages[name] = round(time.time() - t, 1)
        log(f"stage {name} done in {stages[name]} s")

    t = time.time()
    carbon_version = a.gca
    if a.carbon_metadata:
        meta = json.load(open(a.carbon_metadata))
        carbon_version = (f"{a.gca} (Carbon-A post-processing {meta.get('postprocessing', {}).get('revision', '?')}, "
                          f"model {meta.get('probability_model', {}).get('repository', '?')})")
    dbs = {"carbon_gff": _db("Carbon-A annotations.gff3", a.carbon_gff, carbon_version),
           "carbon_proteins": _db("Carbon-A proteins.faa", a.carbon_proteins, carbon_version)}
    loci = io.read_carbon_gff(a.carbon_gff)
    reference, gb2rs = None, None
    if a.ref_gff:
        reference = io.read_reference_gff(a.ref_gff)
        gb2rs = io.read_assembly_report(a.ref_report)
        dbs["ref_gff"] = _db("RefSeq annotation GFF3", a.ref_gff, a.ref_name)
        dbs["ref_proteins"] = _db("RefSeq proteome", a.ref_proteins, a.ref_name)
        dbs["ref_report"] = _db("NCBI assembly report", a.ref_report, a.ref_name)
    cands, counts = select(loci, reference, gb2rs)
    proteins = io.read_fasta(a.carbon_proteins, keep={c.locus.id for c in cands})
    log(f"read inputs: {counts}; proteins found {len(proteins)}")
    faa = os.path.join(raw, "candidates.faa")
    with open(faa + ".tmp", "w") as f:
        for c in cands:
            if c.locus.id in proteins:
                f.write(f">{c.locus.id}\n{proteins[c.locus.id].rstrip('*')}\n")
    os.replace(faa + ".tmp", faa)
    stage("inputs_and_candidates", t)

    tool_versions = {"gene-evidence": __version__}
    tv_path = os.path.join(raw, "tools.json")
    if os.path.exists(tv_path):
        tool_versions.update(json.load(open(tv_path)))
    dlog = open(os.path.join(a.out, "tools.log"), "a")

    def save_versions():
        with open(tv_path, "w") as f:
            json.dump(tool_versions, f, indent=1, sort_keys=True)

    sprot_tsv = os.path.join(raw, "sprot_hits.tsv")
    if a.swissprot:
        dbs["swissprot"] = _db("UniProtKB/Swiss-Prot", a.swissprot, a.swissprot_release)
        if not os.path.exists(sprot_tsv):
            t = time.time()
            tool_versions["diamond"] = tools.diamond_version(a.diamond)
            save_versions()
            log("DIAMOND vs Swiss-Prot: start")
            tools.run_diamond(a.diamond, faa, a.swissprot, os.path.join(raw, "sprot"), sprot_tsv,
                              tools.SWISSPROT_PARAMS, a.threads, dlog)
            stage("diamond_swissprot", t)
    para_tsv = os.path.join(raw, "paralog_hits.tsv")
    if reference is not None and not os.path.exists(para_tsv):
        t = time.time()
        tool_versions["diamond"] = tools.diamond_version(a.diamond)
        save_versions()
        log("DIAMOND vs reference proteome: start")
        tools.run_diamond(a.diamond, faa, a.ref_proteins, os.path.join(raw, "refprot"), para_tsv,
                          tools.PARALOG_PARAMS, a.threads, dlog)
        stage("diamond_paralogs", t)
    pfam_tsv = os.path.join(raw, "pfam_hits.tsv")
    if a.pfam:
        dbs["pfam"] = _db("Pfam-A HMM library", a.pfam, a.pfam_release)
        if not os.path.exists(pfam_tsv):
            t = time.time()
            tool_versions["pyhmmer"] = tools.pyhmmer_version()
            save_versions()
            log("pyhmmer hmmsearch vs Pfam-A: start")
            tools.run_hmmsearch(a.pfam, faa, pfam_tsv, a.threads)
            stage("hmmsearch_pfam", t)

    t = time.time()
    rna_tx = None
    if a.rna_gtf:
        dbs["rna_gtf"] = _db("user-supplied transcript GTF", a.rna_gtf, "user-supplied")
        rna_tx = rna.read_gtf(a.rna_gtf)
    g = graph.build(cands, counts, reference, proteins,
                    tools.parse_diamond(sprot_tsv) if os.path.exists(sprot_tsv) else None,
                    tools.parse_diamond(para_tsv) if reference is not None else None,
                    tools.parse_pfam(pfam_tsv) if os.path.exists(pfam_tsv) else None,
                    rna_tx, dbs, tool_versions, a.top)
    text = report.render(g, a.gca)
    with open(os.path.join(a.out, "graph.json"), "w") as f:
        f.write(graph.dumps(g))
    with open(os.path.join(a.out, "report.md"), "w") as f:
        f.write(text)
    problems = check.check(text, [n["id"] for n in g["nodes"]])
    stage("graph_and_report", t)
    with open(os.path.join(a.out, "versions.json"), "w") as f:
        json.dump({"tools": tool_versions, "inputs_and_databases": dbs,
                   "python": platform.python_version(), "parameters": {
                       "swissprot": tools.SWISSPROT_PARAMS, "paralogs": tools.PARALOG_PARAMS,
                       "pfam": tools.PFAM_PARAMS, "threads": a.threads, "top": a.top}},
                  f, indent=1, sort_keys=True)
        f.write("\n")
    timings = {"wall_seconds": round(time.time() - t0, 1), "stages_seconds": stages,
               "peak_rss_mb_self": _rss_mb(resource.RUSAGE_SELF),
               "peak_rss_mb_largest_child": _rss_mb(resource.RUSAGE_CHILDREN),
               "threads": a.threads, "cpu_count": os.cpu_count(), "platform": platform.platform(),
               "note": "a stage missing from stages_seconds was resumed from an existing raw table"}
    with open(os.path.join(a.out, "timings.json"), "w") as f:
        json.dump(timings, f, indent=1, sort_keys=True)
        f.write("\n")
    log(f"done: {counts['candidates']} candidates, {len(g['ranking'])} in report; "
        f"report check {'PASS' if not problems else 'FAIL'}; {timings}")
    for p in problems[:20]:
        log(p)
    return 1 if problems else 0


def check_report(a):
    g = json.load(open(a.graph))
    problems = check.check(open(a.report).read(), [n["id"] for n in g["nodes"]])
    for p in problems:
        print(p)
    print("PASS" if not problems else f"FAIL: {len(problems)} problem(s)")
    return 1 if problems else 0


def do_fetch(a):
    if a.what == "carbon":
        res = fetch.carbon(a.gca, a.division, a.out)
    elif a.what == "refseq":
        res = fetch.refseq_release(a.release_url, a.gcf, a.out)
    elif a.what == "swissprot":
        res = fetch.swissprot(a.out)
    else:
        res = fetch.pfam(a.out, a.parts)
    print(json.dumps(res, indent=1, sort_keys=True))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(prog="gene-evidence", description="gene-evidence: evidence for Carbon-A "
                                "loci the reference annotation does not have.")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    an = sub.add_parser("analyze", help="candidates, evidence graph and report")
    an.add_argument("--gca", required=True, help="GenBank assembly accession (recorded)")
    an.add_argument("--carbon-gff", required=True)
    an.add_argument("--carbon-proteins", required=True)
    an.add_argument("--carbon-metadata", help="Carbon-A metadata.json (optional, recorded)")
    an.add_argument("--ref-gff", help="RefSeq *_genomic.gff.gz; omit for no reference")
    an.add_argument("--ref-proteins")
    an.add_argument("--ref-report", help="RefSeq *_assembly_report.txt")
    an.add_argument("--ref-name", default="", help="reference release name (recorded)")
    an.add_argument("--swissprot", help="uniprot_sprot.fasta(.gz)")
    an.add_argument("--swissprot-release", default="")
    an.add_argument("--pfam", help="Pfam-A.hmm")
    an.add_argument("--pfam-release", default="")
    an.add_argument("--rna-gtf", help="optional transcript GTF")
    an.add_argument("--diamond", default="diamond")
    an.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    an.add_argument("--top", type=int, default=100, help="candidates kept in graph and report")
    an.add_argument("--out", required=True)
    an.set_defaults(fn=analyze)

    ck = sub.add_parser("check-report", help="fail if a report sentence lacks a valid graph reference")
    ck.add_argument("report")
    ck.add_argument("graph")
    ck.set_defaults(fn=check_report)

    fe = sub.add_parser("fetch", help="download and verify inputs")
    fe.add_argument("what", choices=["carbon", "refseq", "swissprot", "pfam"])
    fe.add_argument("--out", required=True)
    fe.add_argument("--gca")
    fe.add_argument("--division", help="Carbon-A division, e.g. vertebrate_mammalian")
    fe.add_argument("--release-url", help="NCBI annotation-release directory URL")
    fe.add_argument("--gcf")
    fe.add_argument("--parts", type=int, default=16)
    fe.set_defaults(fn=do_fetch)

    a = p.parse_args(argv)
    if a.cmd == "analyze" and a.ref_gff and not (a.ref_proteins and a.ref_report):
        p.error("--ref-gff needs --ref-proteins and --ref-report")
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
