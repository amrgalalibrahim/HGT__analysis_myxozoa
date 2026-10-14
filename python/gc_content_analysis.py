#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
	Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and Evolutionary Implications		
		Amr Galal Abdelraheem Ibrahim; Edson Aparecido Adriano	

	GC content analysis for HGT candidates and ±flank regions.

Implements manuscript-described analysis:
- Compute GC% of candidate interval
- Compute GC% of left flank, right flank, and combined flank (±N bp)
- Report delta GC% (candidate - combined flank)
- Flag deviations that may indicate recent HGT before amelioration

Inputs:
- Genome assembly FASTA
- Candidate coordinates TSV (BED-like)

Output:
- TSV with GC metrics per candidate

Python: 3.8+

Creator: Ibrahim A. G. A. (amrgalalibrahim@gmail.com)
"""

import argparse
import csv
import math
from typing import Dict, List, Optional, Tuple

from utils_fasta import read_fasta, gc_fraction, slice_safe


def parse_candidates(path: str):
    """Parse BED-like TSV candidates.

    Required columns:
      - candidate_id
      - scaffold
      - start (0-based inclusive)
      - end   (0-based exclusive)
    Optional:
      - strand
    """
    with open(path, "r", encoding="utf-8") as fh:
        # Try header-aware parsing
        first = fh.readline()
        if not first:
            raise ValueError("Empty candidates file")
        fh.seek(0)
        dialect = csv.excel_tab
        reader = csv.DictReader(fh, dialect=dialect)
        required = {"candidate_id", "scaffold", "start", "end"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "candidates.tsv must be tab-delimited with header including: "
                "candidate_id, scaffold, start, end (optional strand)"
            )
        for row in reader:
            cid = row["candidate_id"].strip()
            scf = row["scaffold"].strip()
            start = int(row["start"])
            end = int(row["end"])
            strand = row.get("strand", "").strip() or "."
            yield cid, scf, start, end, strand


def pct(x: float) -> float:
    if math.isnan(x):
        return float("nan")
    return 100.0 * x


def main():
    ap = argparse.ArgumentParser(description="GC content analysis for HGT candidates and flanks.")
    ap.add_argument("--genome", required=True, help="Genome assembly FASTA")
    ap.add_argument("--candidates", required=True, help="Candidates TSV (BED-like with header)")
    ap.add_argument("--flank", type=int, default=10000, help="Flank size in bp (default 10000)")
    ap.add_argument("--flag_delta", type=float, default=5.0,
                    help="Flag if |candidate_gc - flank_gc| >= this many percentage points (default 5.0)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    genome = read_fasta(args.genome)

    out_fields = [
        "candidate_id", "scaffold", "start", "end", "strand",
        "len_candidate", "len_left_flank", "len_right_flank",
        "gc_candidate_pct", "gc_left_flank_pct", "gc_right_flank_pct", "gc_flanks_combined_pct",
        "delta_candidate_minus_flanks_pct", "flag_deviation"
    ]

    with open(args.out, "w", encoding="utf-8", newline="") as oh:
        writer = csv.DictWriter(oh, fieldnames=out_fields, dialect=csv.excel_tab)
        writer.writeheader()

        for cid, scf, start, end, strand in parse_candidates(args.candidates):
            if scf not in genome:
                raise KeyError(f"Scaffold '{scf}' not found in genome FASTA.")

            seq = genome[scf]
            cand = slice_safe(seq, start, end)
            left = slice_safe(seq, start - args.flank, start)
            right = slice_safe(seq, end, end + args.flank)
            flanks = left + right

            gc_c = pct(gc_fraction(cand))
            gc_l = pct(gc_fraction(left))
            gc_r = pct(gc_fraction(right))
            gc_f = pct(gc_fraction(flanks))

            if math.isnan(gc_c) or math.isnan(gc_f):
                delta = float("nan")
                flag = "NA"
            else:
                delta = gc_c - gc_f
                flag = "YES" if abs(delta) >= args.flag_delta else "NO"

            writer.writerow({
                "candidate_id": cid,
                "scaffold": scf,
                "start": start,
                "end": end,
                "strand": strand,
                "len_candidate": len(cand),
                "len_left_flank": len(left),
                "len_right_flank": len(right),
                "gc_candidate_pct": f"{gc_c:.3f}" if not math.isnan(gc_c) else "NA",
                "gc_left_flank_pct": f"{gc_l:.3f}" if not math.isnan(gc_l) else "NA",
                "gc_right_flank_pct": f"{gc_r:.3f}" if not math.isnan(gc_r) else "NA",
                "gc_flanks_combined_pct": f"{gc_f:.3f}" if not math.isnan(gc_f) else "NA",
                "delta_candidate_minus_flanks_pct": f"{delta:.3f}" if not math.isnan(delta) else "NA",
                "flag_deviation": flag
            })


if __name__ == "__main__":
    main()
