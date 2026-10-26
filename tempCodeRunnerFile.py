@'
"""Summarise sequences in a FASTA file: length, GC content, base composition."""

import argparse
import csv
import sys
from collections import Counter


def read_fasta(handle):
    """Yield (header, sequence) for each record in an open FASTA file."""
    header, chunks = None, []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        else:
            if header is None:
                raise ValueError("Not a FASTA file: sequence found before first '>' header")
            chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def summarise(header, seq):
    """Return a dict of basic properties for one sequence."""
    seq = seq.upper().replace("U", "T")
    counts = Counter(seq)
    length = len(seq)
    a, c, g, t, n = (counts[b] for b in "ACGTN")
    other = length - (a + c + g + t + n)
    acgt = a + c + g + t
    gc = 100 * (g + c) / acgt