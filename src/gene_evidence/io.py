"""Readers for the input files: Carbon-A GFF3/FASTA, RefSeq GFF3, NCBI assembly report."""
from __future__ import annotations

import collections
import gzip
import hashlib
from dataclasses import dataclass, field
from urllib.parse import unquote

from .intervals import IntervalIndex, union


def open_text(path):
    with open(path, "rb") as f:
        magic = f.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def attrs(column9):
    out = {}
    for part in column9.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = unquote(v)
    return out


def gff_rows(path):
    with open_text(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) >= 9:
                yield c


def read_fasta(path, keep=None):
    """id -> sequence (first word of the header). With `keep`, only those ids."""
    seqs, cur, buf = {}, None, []
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                sid = line[1:].split()[0]
                cur = sid if keep is None or sid in keep else None
                buf = []
            elif cur is not None:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


# ---- Carbon-A ---------------------------------------------------------------------------------------

@dataclass
class Locus:
    id: str
    seq: str
    strand: str
    start: int
    end: int
    status: str
    confidence: str
    cds: list = field(default_factory=list)


def read_carbon_gff(path):
    """Carbon-A annotations.gff3(.gz) -> {locus id: Locus}. CDS = union of the locus's CDS segments."""
    loci, tx2gene = {}, {}
    for c in gff_rows(path):
        a = attrs(c[8])
        if c[2] == "gene":
            loci[a["ID"]] = Locus(a["ID"], c[0], c[6], int(c[3]), int(c[4]),
                                  a.get("status", ""), a.get("confidence", ""))
        elif c[2] in ("mRNA", "transcript"):
            tx2gene[a["ID"]] = a["Parent"]
        elif c[2] == "CDS":
            parent = a["Parent"].split(",")[0]
            gene = tx2gene.get(parent, parent)
            loci[gene].cds.append((int(c[3]), int(c[4])))
    for locus in loci.values():
        locus.cds = union(locus.cds)
    return loci


# ---- RefSeq reference --------------------------------------------------------------------------------

@dataclass
class RefGene:
    gene_id: str
    name: str
    seq: str
    strand: str
    start: int
    end: int
    biotype: str


@dataclass
class Reference:
    pc_cds: dict          # (seq, strand) -> IntervalIndex of protein-coding CDS (pseudo=true excluded)
    genes: dict           # seq -> [RefGene] for gene/pseudogene features, sorted by (start, end, gene_id)
    proteins: dict        # protein_id -> {"gene_id", "segments"}
    gene_by_id: dict      # gene_id -> RefGene


def _gene_id(a):
    for x in a.get("Dbxref", "").split(","):
        if x.startswith("GeneID:"):
            return x[len("GeneID:"):]
    return a.get("ID", "")


def read_reference_gff(path):
    cds = collections.defaultdict(list)
    genes = collections.defaultdict(list)
    prot_segments = collections.defaultdict(list)
    prot_gene = {}
    for c in gff_rows(path):
        kind = c[2]
        if kind in ("gene", "pseudogene"):
            a = attrs(c[8])
            biotype = a.get("gene_biotype") or ("pseudogene" if kind == "pseudogene" else "unknown")
            genes[c[0]].append(RefGene(_gene_id(a), a.get("Name", a.get("gene", "")), c[0], c[6],
                                       int(c[3]), int(c[4]), biotype))
        elif kind == "CDS":
            a = attrs(c[8])
            if a.get("pseudo") == "true":
                continue
            cds[(c[0], c[6])].append((int(c[3]), int(c[4])))
            pid = a.get("protein_id")
            if pid:
                prot_segments[pid].append((int(c[3]), int(c[4])))
                prot_gene[pid] = _gene_id(a)
    for seq in genes:
        genes[seq].sort(key=lambda g: (g.start, g.end, g.gene_id))
    gene_by_id = {g.gene_id: g for gs in genes.values() for g in gs}
    proteins = {p: {"gene_id": prot_gene[p], "segments": len(union(s))} for p, s in prot_segments.items()}
    return Reference({k: IntervalIndex(v) for k, v in cds.items()}, dict(genes), proteins, gene_by_id)


def read_assembly_report(path):
    """GenBank accession -> RefSeq accession, for rows whose Relationship is '=' (identical sequence)."""
    out = {}
    with open_text(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            # Sequence-Name Role Assigned-Molecule Location/Type GenBank-Accn Relationship RefSeq-Accn ...
            if len(c) > 6 and c[5] == "=" and c[4] != "na" and c[6] != "na":
                out[c[4]] = c[6]
    return out
