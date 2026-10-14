#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
	Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and Evolutionary Implications		
		Amr Galal Abdelraheem Ibrahim; Edson Aparecido Adriano	

		FASTA utilities (pure Python, no external dependencies).

Designed for Python 3.8 compatibility.

Creator: Ibrahim A. G. A. (amrgalalibrahim@gmail.com)
"""

from typing import Dict, Iterator, Tuple, TextIO, Optional


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA into a dict: {record_id: sequence}.

    Notes:
      - Record id is taken as the first token after '>' (up to first whitespace).
      - Sequences are uppercased and concatenated across lines.
    """
    seqs: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as fh:
        rec_id: Optional[str] = None
        chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rec_id is not None:
                    seqs[rec_id] = "".join(chunks).upper()
                rec_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if rec_id is not None:
            seqs[rec_id] = "".join(chunks).upper()
    return seqs


def gc_fraction(seq: str) -> float:
    """Compute GC fraction among A/C/G/T only. Returns NaN if no A/C/G/T."""
    seq = seq.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    denom = a + c + g + t
    if denom == 0:
        return float("nan")
    return (c + g) / denom


def slice_safe(seq: str, start: int, end: int) -> str:
    """Safe slice (clamps to [0, len(seq)])."""
    n = len(seq)
    s = max(0, min(n, start))
    e = max(0, min(n, end))
    if e < s:
        s, e = e, s
    return seq[s:e]
