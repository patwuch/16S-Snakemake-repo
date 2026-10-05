#!/usr/bin/env python3
"""The factorial report's abundance table, computed from the ALDEx3 tables.

One row per report comparison: taxa with q < 0.05 in the primary fit
(prevalence >= 6, gamma = 0.5), the largest increase and decrease among them,
their smallest q, and the same count under the stricter filter (>= 12) and at
gamma = 1. Shared by the report text and check_report_factorial.py.

Usage: python3 aldex3_rows_factorial.py   (prints the rows)
"""
import os

import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(THIS_DIR), "comparison-Treatment-x-Radiation-factorial", "tables")
ALPHA = 0.05

# (question, report label, ALDEx3 row, stratum)
ROWS = [
    ("Do the treatments differ in unirradiated mice?", "IFA vs Control", "IFA vs Control", "unirradiated"),
    ("Do the treatments differ in unirradiated mice?", "Inosine vs Control", "Inosine vs Control", "unirradiated"),
    ("Do the treatments differ in unirradiated mice?", "IFA vs Inosine", "IFA vs Inosine", "unirradiated"),
    ("Do the treatments differ in irradiated mice?", "IFA+IR vs IR", "IFA vs Control", "irradiated"),
    ("Do the treatments differ in irradiated mice?", "Inosine+IR vs IR", "Inosine vs Control", "irradiated"),
    ("Do the treatments differ in irradiated mice?", "IFA+IR vs Inosine+IR", "IFA vs Inosine", "irradiated"),
    ("Does radiation change the microbiome within each treatment?", "Control", "Radiation within Control", ""),
    ("Does radiation change the microbiome within each treatment?", "IFA", "Radiation within IFA", ""),
    ("Does radiation change the microbiome within each treatment?", "Inosine", "Radiation within Inosine", ""),
    ("Does radiation's effect differ between treatments?", "IFA vs Control", "Interaction: IFA vs Control", ""),
    ("Does radiation's effect differ between treatments?", "Inosine vs Control", "Interaction: Inosine vs Control", ""),
    ("Does radiation's effect differ between treatments?", "IFA vs Inosine", "Interaction: IFA vs Inosine", ""),
]


def table(minpres=6, gamma="0.5", extra=""):
    d = pd.read_csv(os.path.join(TAB, f"ALDEx3_factorial_Greengenes2_Day17_minpres{minpres}_gamma{gamma}{extra}.tsv"),
                    sep="\t", keep_default_na=False)
    for c in ("log2FC", "p", "q", "ci_lo", "ci_hi"):
        d[c] = d[c].astype(float)
    return d


def pick(d, row, stratum):
    return d[(d.row == row) & (d.stratum == stratum)]


def rows():
    main, strict, wide = table(), table(12), table(6, "1")
    out = []
    for question, lab, row, st in ROWS:
        d = pick(main, row, st)
        h = d[d.q < ALPHA]
        up, down = h[h.log2FC > 0], h[h.log2FC < 0]
        out.append(dict(
            question=question, label=lab, row=row, stratum=st, n=len(h),
            max_up=up.log2FC.max() if len(up) else None,
            max_down=-down.log2FC.min() if len(down) else None,
            min_q=h.q.min() if len(h) else None,
            n_strict=int((pick(strict, row, st).q < ALPHA).sum()),
            n_gamma1=int((pick(wide, row, st).q < ALPHA).sum()),
            nominal=int((d.p < ALPHA).sum()), tested=len(d),
            hits=h.sort_values("q")[["taxon", "log2FC", "q"]].to_dict("records"),
            all=d))
    return out


if __name__ == "__main__":
    for r in rows():
        print(f"{r['question'][:30]:30} {r['label']:22} n={r['n']} up={r['max_up']} down={r['max_down']} "
              f"minq={r['min_q']} strict={r['n_strict']} g1={r['n_gamma1']} nominal={r['nominal']}/{r['tested']}")
        for h in r["hits"]:
            print(f"      {h['taxon'].split(';')[-2:]} {h['log2FC']:+.2f} q={h['q']:.3g}")
