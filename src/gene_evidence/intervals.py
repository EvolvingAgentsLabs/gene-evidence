"""Closed integer intervals (GFF coordinates: 1-based, inclusive)."""
from __future__ import annotations

import bisect


def union(intervals):
    """Sorted, merged copy of `intervals` (touching intervals are merged)."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


class IntervalIndex:
    """Sorted disjoint intervals with fast overlap queries."""

    def __init__(self, intervals):
        u = union(intervals)
        self.starts = [s for s, _ in u]
        self.ends = [e for _, e in u]

    def overlap_bp(self, query):
        """Number of bp of the (merged) query intervals covered by this index."""
        total = 0
        for s, e in union(query):
            i = bisect.bisect_left(self.ends, s)  # first interval ending at or after s
            while i < len(self.starts) and self.starts[i] <= e:
                total += min(e, self.ends[i]) - max(s, self.starts[i]) + 1
                i += 1
        return total


def span_overlaps(a_start, a_end, b_start, b_end):
    return a_start <= b_end and b_start <= a_end


def distance(a_start, a_end, b_start, b_end):
    """Gap in bp between two intervals; 0 if they overlap or touch."""
    if span_overlaps(a_start, a_end, b_start, b_end):
        return 0
    return max(b_start - a_end, a_start - b_end) - 1
