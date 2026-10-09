"""Candidate selection: Carbon-A loci whose CDS overlaps no reference protein-coding CDS (same strand)."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Candidate:
    locus: object        # io.Locus
    ref_seq: str         # RefSeq sequence name ("" without a reference)


def select(loci, reference=None, gb2rs=None):
    """Returns (candidates sorted by locus id, counts).

    Without a reference every locus is a candidate. With one, loci on a GenBank sequence that has no
    identical RefSeq twin are dropped and counted."""
    counts = {"carbon_loci": len(loci), "dropped_no_identical_twin": 0,
              "overlap_reference_cds": 0, "candidates": 0}
    out = []
    for lid in sorted(loci):
        locus = loci[lid]
        if reference is None:
            out.append(Candidate(locus, ""))
            continue
        rs = gb2rs.get(locus.seq)
        if rs is None:
            counts["dropped_no_identical_twin"] += 1
            continue
        idx = reference.pc_cds.get((rs, locus.strand))
        if idx is not None and idx.overlap_bp(locus.cds) > 0:
            counts["overlap_reference_cds"] += 1
            continue
        out.append(Candidate(locus, rs))
    counts["candidates"] = len(out)
    return out, counts
