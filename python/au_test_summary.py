#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
au_test_summary.py
==================
Summarize the Approximately Unbiased (AU) topology tests performed with
IQ-TREE v2.1.3 (`-au` flag).  For each candidate gene, read the
*.iqtree log and extract the AU p-value for the alternative topology that
constrains the myxozoan clade as sister to all other eukaryotes (not
specifically sister to fungi).  Output a single summary TSV.

Usage:
    python au_test_summary.py --tree_dir <dir> --out <file.tsv>
"""
import argparse
import csv
import re
from pathlib import Path


def parse_iqtree_log(path):
    """Extract the AU p-value for the constrained topology from a *.iqtree log."""
    au_p = None
    constrained_logL = None
    unconstrained_logL = None
    in_au = False
    for line in path.read_text().splitlines():
        if "USER TREES" in line:
            in_au = False
        if "Approximately Unbiased" in line or "au" in line.lower():
            in_au = True
        if in_au:
            m = re.search(r"(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
            if m:
                # Last numeric token is usually the AU p
                tokens = m.groups()
                try:
                    au_p = float(tokens[-1])
                except ValueError:
                    pass
        if "BEST SCORE FOUND" in line:
            try:
                unconstrained_logL = float(line.split(":")[-1].strip())
            except ValueError:
                pass
    return au_p, unconstrained_logL, constrained_logL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree_dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    tree_dir = Path(args.tree_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, dialect=csv.excel_tab)
        w.writerow([
            "candidate", "AU_pvalue_alternative_topology",
            "logL_unconstrained", "logL_constrained", "reject_alt_p_lt_0.05"
        ])
        for log in sorted(tree_dir.glob("*ML.iqtree")):
            cand = log.name.replace("_ML.iqtree", "")
            au_p, l_un, l_co = parse_iqtree_log(log)
            reject = "YES" if (au_p is not None and au_p < 0.05) else "NO"
            w.writerow([cand,
                        f"{au_p:.4f}" if au_p is not None else "NA",
                        f"{l_un:.2f}" if l_un is not None else "NA",
                        f"{l_co:.2f}" if l_co is not None else "NA",
                        reject])
    print(f"Wrote AU-test summary to {out_path}")


if __name__ == "__main__":
    main()
