"""External tools: DIAMOND (blastp) and pyhmmer (hmmsearch). Runners write raw tables; parsers read them.

Each runner writes to a temporary file and renames it when done, so a table that exists is complete.
`gene-evidence analyze` skips a stage whose table already exists, which makes an interrupted run resumable."""
from __future__ import annotations

import os
import re
import subprocess

DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore",
                  "qcovhsp", "stitle"]
SWISSPROT_PARAMS = {"mode": "blastp, default sensitivity", "evalue": 1e-5, "max_target_seqs": 25}
PARALOG_PARAMS = {"mode": "blastp, default sensitivity", "evalue": 1e-5, "max_target_seqs": 25}
PFAM_PARAMS = {"program": "hmmsearch", "bit_cutoffs": "gathering", "reported": "included domains"}


def diamond_version(diamond):
    out = subprocess.run([diamond, "version"], capture_output=True, text=True, check=True).stdout
    m = re.search(r"(\d+\.\d+\.\d+)", out)
    return m.group(1) if m else out.strip()


def run_diamond(diamond, query_faa, db_fasta, db_path, out_tsv, params, threads, log=None):
    if not os.path.exists(db_path + ".dmnd"):
        subprocess.run([diamond, "makedb", "--in", db_fasta, "-d", db_path, "-p", str(threads)],
                       check=True, stdout=log, stderr=log)
    tmp = out_tsv + ".tmp"
    subprocess.run([diamond, "blastp", "-q", query_faa, "-d", db_path, "-o", tmp,
                    "--evalue", str(params["evalue"]), "-k", str(params["max_target_seqs"]),
                    "-p", str(threads), "--outfmt", "6", *DIAMOND_FIELDS],
                   check=True, stdout=log, stderr=log)
    os.replace(tmp, out_tsv)


def parse_diamond(path):
    """qseqid -> best hit (highest bitscore; ties by sseqid), as a dict of typed fields."""
    best = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < len(DIAMOND_FIELDS):
                continue
            h = dict(zip(DIAMOND_FIELDS, c))
            for k in ("pident", "evalue", "bitscore", "qcovhsp"):
                h[k] = float(h[k])
            for k in ("length", "qlen", "slen"):
                h[k] = int(h[k])
            cur = best.get(h["qseqid"])
            if cur is None or (-h["bitscore"], h["sseqid"]) < (-cur["bitscore"], cur["sseqid"]):
                best[h["qseqid"]] = h
    return best


def swissprot_title(stitle):
    """'sp|P1|X_HUMAN Name OS=Homo sapiens OX=9606 GN=..' -> (name, organism, taxid)."""
    name = re.sub(r"^\S+\s*", "", stitle)
    name = name.split(" OS=")[0].strip()
    os_m = re.search(r" OS=(.+?)(?= [A-Z]{2}=|$)", stitle)
    ox_m = re.search(r" OX=(\d+)", stitle)
    return name, os_m.group(1) if os_m else "", ox_m.group(1) if ox_m else ""


def pyhmmer_version():
    import pyhmmer
    return pyhmmer.__version__


def _text(x):
    return x.decode() if isinstance(x, bytes) else (x or "")


def run_hmmsearch(pfam_hmm, query_faa, out_tsv, cpus):
    """Pfam-A HMMs (queries) against the candidate proteins (targets), gathering thresholds."""
    import pyhmmer
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(query_faa, digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()
    rows = []
    with pyhmmer.plan7.HMMFile(pfam_hmm) as hf:
        for hits in pyhmmer.hmmsearch(hf, seqs, cpus=cpus, bit_cutoffs="gathering"):
            q = hits.query
            acc, name = _text(getattr(q, "accession", "")), _text(getattr(q, "name", ""))
            for hit in hits.included:
                for d in hit.domains.included:
                    rows.append((_text(hit.name), acc, name, f"{d.score:.1f}", f"{d.i_evalue:.2g}",
                                 str(d.env_from), str(d.env_to)))
    rows.sort(key=lambda r: (r[0], int(r[5]), r[1]))
    tmp = out_tsv + ".tmp"
    with open(tmp, "w") as f:
        f.write("locus\tpfam_acc\tpfam_name\tscore\ti_evalue\tenv_from\tenv_to\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    os.replace(tmp, out_tsv)


def parse_pfam(path):
    """locus -> list of domain dicts, in table order."""
    out = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            out.setdefault(r["locus"], []).append({
                "acc": r["pfam_acc"], "name": r["pfam_name"], "score": float(r["score"]),
                "i_evalue": float(r["i_evalue"]), "env_from": int(r["env_from"]),
                "env_to": int(r["env_to"])})
    return out
