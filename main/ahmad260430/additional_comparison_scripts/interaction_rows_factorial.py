#!/usr/bin/env python3
"""The interaction rows of comparison-Treatment-x-Radiation-factorial/report.html.

One row per treatment pair, from plot_interaction_forest.py's panels and joint
q (one Benjamini-Hochberg pass over every taxon in all three panels). A hit is
q < 0.05 and ANCOM-BC2's passed_ss. check_report_factorial.py compares the
report against these values; run this directly to print them.

Usage, inside the qiime2 container:
    python3 interaction_rows_factorial.py
"""
import math
import os

import pandas as pd

from plot_interaction_forest import joint_q, load_panels

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TAB_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Treatment-x-Radiation-factorial", "tables")
LABEL = {"TreatmentIFA:RadiationIR": "IFA vs Control", "TreatmentInosine:RadiationIR": "Inosine vs Control",
         "IFA vs Inosine": "IFA vs Inosine"}


def short(t):
    parts = [p.strip() for p in t.split(";") if p.strip() and not p.strip().endswith("__")]
    return parts[-1] if parts else t


def rows(suffix="_Day17_minpres6"):
    tables, panels = load_panels(TAB_DIR, f"ANCOMBC2_factorial_Greengenes2{suffix}")
    qj = joint_q(tables, panels)
    for panel, (fit, coef, sign) in panels.items():
        df = tables[fit].set_index("taxon")
        hits = [(t, sign * df.loc[t, f"lfc_{coef}"] / math.log(2), q) for t, q in qj[panel].items()
                if q < 0.05 and str(df.loc[t, f"passed_ss_{coef}"]).upper() == "TRUE"]
        up = [h for h in hits if h[1] > 0]
        down = [h for h in hits if h[1] < 0]
        nominal = sum(1 for t in qj[panel] if df.loc[t, f"p_{coef}"] < 0.05) - len(hits)
        yield dict(label=LABEL[panel], n=len(hits), n_tested=len(qj[panel]), nominal=nominal,
                   max_up=max((h[1] for h in up), default=None), max_down=max((-h[1] for h in down), default=None),
                   min_q=min((h[2] for h in hits), default=None),
                   hits=sorted((short(t), round(l, 2), q) for t, l, q in hits))


if __name__ == "__main__":
    for suffix in ("_Day17_minpres6", "_Day17_minpres12"):
        for r in rows(suffix):
            print(suffix, r["label"], f"n={r['n']} of {r['n_tested']} (nominal-only {r['nominal']})", r["hits"])
