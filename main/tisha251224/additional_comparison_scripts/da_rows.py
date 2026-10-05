#!/usr/bin/env python3
"""The MaAsLin3 differential-abundance rows report.html shows, from run_maaslin3_rows.py's tables.

One row per report row (question x contrast): features tested, calls at
q < 0.05, calls robust to dropping any one mouse, largest change in share in
each direction among calls, smallest q, the typical (median) taxon's change,
and the status chip. Also each call's details and whether it is a call at the
>= 6-sample filter. check_report.py and gen_da_html.py read these; run this
directly to print them.

Usage, inside the qiime2 container:
    python3 da_rows.py
"""
import os
import re

import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
MA_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial", "tables", "maaslin3")

Q_DIET0 = "Do the diets differ in unirradiated mice?"
Q_DIET1 = "Do the diets differ in irradiated mice?"
Q_IR = "Does irradiation change the microbiome within each diet?"
Q_INT = "Does irradiation's effect differ between diets?"
# (question, label, first-named group, second, run_maaslin3_rows row, stratum, .res suffix or None)
ROWS = [
    (Q_DIET0, "Inulin vs Cellulose", "Inulin", "Cellulose", "Inulin vs Cellulose", "Control", "_DietInulin"),
    (Q_DIET0, "RS vs Cellulose", "RS", "Cellulose", "RS vs Cellulose", "Control", "_DietRS"),
    (Q_DIET0, "RS vs Inulin", "RS", "Inulin", "RS vs Inulin", "Control", "_refInulin_DietRS"),
    (Q_DIET1, "Inulin+IR vs Cellulose+IR", "Inulin+IR", "Cellulose+IR", "Inulin vs Cellulose", "IR", "_condrefIR_DietInulin"),
    (Q_DIET1, "RS+IR vs Cellulose+IR", "RS+IR", "Cellulose+IR", "RS vs Cellulose", "IR", "_condrefIR_DietRS"),
    (Q_DIET1, "RS+IR vs Inulin+IR", "RS+IR", "Inulin+IR", "RS vs Inulin", "IR", "_condrefIR_refInulin_DietRS"),
    (Q_IR, "Cellulose+IR vs Cellulose", "Cellulose+IR", "Cellulose", "IR within Cellulose", "", "_ConditionIR"),
    (Q_IR, "Inulin+IR vs Inulin", "Inulin+IR", "Inulin", "IR within Inulin", "", "_refInulin_ConditionIR"),
    (Q_IR, "RS+IR vs RS", "RS+IR", "RS", "IR within RS", "", "_refRS_ConditionIR"),
    (Q_INT, "Inulin vs Cellulose", "Inulin", "Cellulose", "Interaction: Inulin vs Cellulose", "", None),
    (Q_INT, "RS vs Cellulose", "RS", "Cellulose", "Interaction: RS vs Cellulose", "", None),
    (Q_INT, "RS vs Inulin", "RS", "Inulin", "Interaction: RS vs Inulin", "", None),
]
N_LOO = 19


def load(min_samples):
    d = pd.read_csv(os.path.join(MA_DIR, f"MaAsLin3_factorial_Greengenes2_min{min_samples}_rows.tsv"), sep="\t")
    d["stratum"] = d.stratum.fillna("")
    return d[d.model == "abundance"]


def species_named(taxon):
    last = [p for p in taxon.split(";") if not p.endswith("__")][-1]
    return bool(re.match(r"s__.", last))


def rows():
    a, b = load(2), load(6)
    b6 = b.set_index(["row", "stratum", "feature"])
    out = []
    for question, label, first, second, row, stratum, res in ROWS:
        d = a[(a.row == row) & (a.stratum == stratum)]
        fitted = d[d.p.notna()]
        calls = fitted[fitted.call].copy()
        calls["call_at_6"] = [bool(b6.loc[(row, stratum, f), "call"]) if (row, stratum, f) in b6.index else False
                              for f in calls.feature]
        n_rob = int(calls.robust.sum())
        chip = ("robust", "ok") if n_rob else (("not robust", "warn") if len(calls) else ("below threshold", "warn"))
        up, down = calls.coef[calls.coef > 0], calls.coef[calls.coef < 0]
        out.append(dict(
            question=question, label=label, first=first, second=second, res=res,
            fitted=len(fitted), n=len(calls), n_robust=n_rob,
            max_up=up.max() if len(up) else None, max_down=-down.min() if len(down) else None,
            min_q=calls.q.min() if len(calls) else None, median=fitted.coef.median(),
            chip=chip, calls=calls.sort_values("coef", ascending=False),
            n_at_6=len(b[(b.row == row) & (b.stratum == stratum) & b.call])))
    return out


if __name__ == "__main__":
    for r in rows():
        print(f"{r['question'][:28]:28s} {r['label']:26s} tested {r['fitted']:3d} calls {r['n']:2d} robust {r['n_robust']} "
              f"up {r['max_up']} down {r['max_down']} minq {r['min_q']} median {r['median']:+.2f} {r['chip'][0]} (>=6: {r['n_at_6']})")
