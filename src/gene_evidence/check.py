"""Report-claim checker (acceptance criterion (b)).

A claim unit is a table data row, or a sentence of any other non-heading line. Each unit must carry at
least one [node-id] reference, and every reference must resolve to a node of the graph."""
from __future__ import annotations

import re

REF = re.compile(r"\[([a-z_]+:[^\[\]\s]+)\]")
LABEL = re.compile(r"^\*\*[^*]+\*\*$")
SEPARATOR = re.compile(r"^\|[\s:|-]+\|$")
SENTENCE_END = re.compile(r"(?<=[.!?])\s+")


def claim_units(text):
    lines = text.splitlines()
    for k, raw in enumerate(lines):
        line = raw.strip()
        if not line or line.startswith("#") or LABEL.match(line) or SEPARATOR.match(line):
            continue
        if line.startswith("|"):
            nxt = lines[k + 1].strip() if k + 1 < len(lines) else ""
            if SEPARATOR.match(nxt):  # table header
                continue
            yield k + 1, line
            continue
        if line.startswith("- "):
            line = line[2:]
        for sentence in SENTENCE_END.split(line):
            if sentence.strip():
                yield k + 1, sentence.strip()


def check(report_text, node_ids):
    """-> list of problems (empty when the report passes)."""
    node_ids = set(node_ids)
    problems = []
    for lineno, unit in claim_units(report_text):
        refs = REF.findall(unit)
        if not refs:
            problems.append(f"line {lineno}: uncited claim: {unit[:120]}")
        for r in refs:
            if r not in node_ids:
                problems.append(f"line {lineno}: reference [{r}] is not a graph node")
    return problems
