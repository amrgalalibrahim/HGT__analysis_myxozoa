#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
salmon_to_counts.py
===================
Aggregate Salmon `quant.sf` files into a single counts matrix
(genes x samples) suitable for DESeq2 input.

Usage:
    python salmon_to_counts.py \\
        --quant_dir <dir_with_quant_subdirs> \\
        --samples   samples.tsv \\
        --out       counts.tsv

`--samples` is a TSV with at least the column `sample`.
The script looks for `<quant_dir>/<sample>/quant.sf`.
"""
import argparse
import csv
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant_dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    quant_dir = Path(args.quant_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.samples) as fh:
        samples = [r["sample"].strip() for r in csv.DictReader(fh, dialect=csv.excel_tab)]

    # Load each sample's quant.sf
    sample_counts = {}
    gene_ids = None
    for s in samples:
        qfile = quant_dir / s / "quant.sf"
        if not qfile.exists():
            print(f"WARNING: {qfile} not found, skipping")
            continue
        per_gene = {}
        with open(qfile) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            name_i = header.index("Name")
            num_i = header.index("NumReads")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                per_gene[parts[name_i]] = int(float(parts[num_i]))
        sample_counts[s] = per_gene
        if gene_ids is None:
            gene_ids = sorted(per_gene.keys())
        else:
            gene_ids = sorted(set(gene_ids) | set(per_gene.keys()))

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, dialect=csv.excel_tab)
        w.writerow(["gene_id"] + samples)
        for g in gene_ids:
            w.writerow([g] + [sample_counts.get(s, {}).get(g, 0) for s in samples])
    print(f"Wrote counts matrix ({len(gene_ids)} genes x {len(samples)} samples) to {out_path}")


if __name__ == "__main__":
    main()
