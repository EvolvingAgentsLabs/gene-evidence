"""Optional RNA evidence from a user-supplied transcript GTF. Absent input -> a recorded no-op.

Stage 1 records it but does not use it for tiering (BRIEF.md, evidence 6). BAM input is out of scope."""
from __future__ import annotations

import collections
import re

from .intervals import span_overlaps
from .io import open_text


def read_gtf(path):
    """(seq, strand) -> {transcript_id: sorted exon list}."""
    tx = collections.defaultdict(lambda: collections.defaultdict(list))
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', c[8])
            if m:
                tx[(c[0], c[6])][m.group(1)].append((int(c[3]), int(c[4])))
    return {k: {t: sorted(v) for t, v in d.items()} for k, d in tx.items()}


def introns(exons):
    return {(a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:])}


def evidence(transcripts, seq_names, strand, start, end, cds):
    """transcripts: read_gtf output or None. seq_names: names the candidate's sequence may carry."""
    if transcripts is None:
        return {"supplied": False}
    cand_introns = introns(sorted(cds))
    over, matched = [], set()
    for seq in seq_names:
        for tid, exons in sorted(transcripts.get((seq, strand), {}).items()):
            if span_overlaps(start, end, exons[0][0], exons[-1][1]):
                over.append(tid)
                matched |= cand_introns & introns(exons)
    return {"supplied": True, "overlapping_transcripts": over,
            "predicted_junctions": len(cand_introns),
            "junctions_matched": sorted([list(j) for j in matched])}
