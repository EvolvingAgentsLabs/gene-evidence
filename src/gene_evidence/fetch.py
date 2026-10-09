"""Downloads with verification. Carbon-A files are checked against the sha256 in Carbon-A's own
metadata.json; NCBI files against the release's md5checksums.txt; Pfam is fetched with parallel HTTP range
requests (EBI throttles single connections) and checked by exact byte count and a full gzip read."""
from __future__ import annotations

import concurrent.futures as cf
import gzip
import hashlib
import json
import os
import re
import shutil
import urllib.request

UA = "EvolvingAgentsLabs research (gene-evidence)"
CARBON = "https://huggingface.co/buckets/HuggingFaceBio/genbank-annotations/resolve/cds/{div}/{gca}/{name}"
SPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
PFAM = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"


def _open(url, headers=None, timeout=300):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA, **(headers or {})}),
                                  timeout=timeout)


def get(url, path, retries=5):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    for attempt in range(retries):
        try:
            with _open(url) as r, open(path + ".part", "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            os.replace(path + ".part", path)
            return path
        except OSError:
            if attempt == retries - 1:
                raise
    return path


def get_ranged(url, path, parts=16):
    """Parallel range download; verifies the exact byte count."""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    with _open(url, timeout=60) as r:
        size = int(r.headers["Content-Length"])
    step = -(-size // parts)

    def one(k):
        a, b = k * step, min(size, (k + 1) * step) - 1
        part = f"{path}.{k:03d}"
        for _ in range(8):
            have = os.path.getsize(part) if os.path.exists(part) else 0
            if have == b - a + 1:
                return
            try:
                with _open(url, {"Range": f"bytes={a + have}-{b}"}) as r, open(part, "ab") as f:
                    shutil.copyfileobj(r, f, 1 << 20)
            except OSError:
                continue
        if os.path.getsize(part) != b - a + 1:
            raise OSError(f"range {k} incomplete")

    with cf.ThreadPoolExecutor(parts) as ex:
        list(ex.map(one, range(parts)))
    with open(path + ".part", "wb") as f:
        for k in range(parts):
            with open(f"{path}.{k:03d}", "rb") as p:
                shutil.copyfileobj(p, f, 1 << 20)
            os.remove(f"{path}.{k:03d}")
    if os.path.getsize(path + ".part") != size:
        raise OSError("size mismatch after range download")
    os.replace(path + ".part", path)
    return path


def _hash(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def gzip_ok(path):
    with gzip.open(path, "rb") as f:
        while f.read(1 << 24):
            pass
    return True


def carbon(gca, division, out):
    os.makedirs(out, exist_ok=True)
    meta = get(CARBON.format(div=division, gca=gca, name="metadata.json"), os.path.join(out, "metadata.json"))
    want = {f["name"]: f["sha256"] for f in json.load(open(meta))["files"]}
    res = {}
    for name in ("annotations.gff3.gz", "proteins.faa.gz"):
        p = get(CARBON.format(div=division, gca=gca, name=name), os.path.join(out, name))
        got = _hash(p, "sha256")
        if got != want[name]:
            raise ValueError(f"{name}: sha256 {got} != metadata {want[name]}")
        res[name] = got
    return res


def refseq_release(release_url, gcf, out):
    """release_url: an NCBI annotation-release directory, e.g.
    https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10116/GCF_036323735.1-RS_2024_02/"""
    os.makedirs(out, exist_ok=True)
    release_url = release_url.rstrip("/") + "/"
    with _open(release_url) as r:
        names = re.findall(r'href="([^"/?][^"]*)"', r.read().decode())
    md5 = {}
    with _open(release_url + "md5checksums.txt") as r:
        for line in r.read().decode().splitlines():
            h, n = line.split(None, 1)
            md5[n.strip().lstrip("./")] = h
    res = {}
    for suffix in ("_genomic.gff.gz", "_protein.faa.gz", "_assembly_report.txt"):
        name = next(n for n in names if n.startswith(gcf + "_") and n.endswith(suffix))
        p = get(release_url + name, os.path.join(out, name))
        if name in md5 and _hash(p, "md5") != md5[name]:
            raise ValueError(f"{name}: md5 mismatch")
        res[name] = {"md5_verified": name in md5, "sha256": _hash(p, "sha256")}
    return res


def swissprot(out):
    os.makedirs(out, exist_ok=True)
    p = get(SPROT + "uniprot_sprot.fasta.gz", os.path.join(out, "uniprot_sprot.fasta.gz"))
    rel = get(SPROT + "reldate.txt", os.path.join(out, "uniprot_reldate.txt"))
    gzip_ok(p)
    text = open(rel).read()
    m = re.search(r"Swiss-Prot Release (\S+) of (\S+)", text)
    return {"file": p, "release": f"{m.group(1)} of {m.group(2)}" if m else text.strip()}


def pfam(out, parts=16):
    os.makedirs(out, exist_ok=True)
    gz = get_ranged(PFAM + "Pfam-A.hmm.gz", os.path.join(out, "Pfam-A.hmm.gz"), parts)
    gzip_ok(gz)
    ver = get(PFAM + "Pfam.version.gz", os.path.join(out, "Pfam.version.gz"))
    text = gzip.open(ver, "rt").read()
    m = re.search(r"Pfam release\s*:\s*(\S+)", text)
    hmm = os.path.join(out, "Pfam-A.hmm")
    if not os.path.exists(hmm):
        with gzip.open(gz, "rb") as fi, open(hmm + ".part", "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 24)
        os.replace(hmm + ".part", hmm)
    return {"file": hmm, "gz_sha256": _hash(gz, "sha256"), "release": m.group(1) if m else text.strip()}
