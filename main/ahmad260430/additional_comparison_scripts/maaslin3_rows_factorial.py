#!/usr/bin/env python3
"""The factorial report's abundance table, computed from the MaAsLin3 rows table.

One row per report comparison: taxa tested, calls (q < 0.05) and those holding in all 24
leave-one-mouse-out refits, the largest increase and decrease among the calls, their
smallest q, and the typical taxon's change (median over the tested taxa). Shared by the
report text and check_report_factorial.py.

Usage: python3 maaslin3_rows_factorial.py   (prints the rows and the table's HTML)
"""
import os

import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
MA = os.path.join(os.path.dirname(THIS_DIR), "comparison-Treatment-x-Radiation-factorial", "tables", "maaslin3")
TABLE = os.path.join(MA, "MaAsLin3_factorial_Greengenes2_Day17_min3_rows.tsv")
ALPHA = 0.05
N_MICE = 24

# (question, report label (baseline on the left), MaAsLin3 row, stratum, swatch)
ROWS = [
    ("Do the treatments differ in unirradiated mice?", "Control vs IFA", "IFA vs Control", "unirradiated", "#90BFF9"),
    ("Do the treatments differ in unirradiated mice?", "Control vs Inosine", "Inosine vs Control", "unirradiated", "#FF8000"),
    ("Do the treatments differ in unirradiated mice?", "Inosine vs IFA", "IFA vs Inosine", "unirradiated", None),
    ("Do the treatments differ in irradiated mice?", "IR vs IFA+IR", "IFA vs Control", "irradiated", "#0000FF"),
    ("Do the treatments differ in irradiated mice?", "IR vs Inosine+IR", "Inosine vs Control", "irradiated", "#FF0000"),
    ("Do the treatments differ in irradiated mice?", "Inosine+IR vs IFA+IR", "IFA vs Inosine", "irradiated", None),
    ("Does radiation change the microbiome within each treatment?", "Control vs IR", "Radiation within Control", "", "#000"),
    ("Does radiation change the microbiome within each treatment?", "IFA vs IFA+IR", "Radiation within IFA", "", "#0000FF"),
    ("Does radiation change the microbiome within each treatment?", "Inosine vs Inosine+IR", "Radiation within Inosine", "", "#FF0000"),
    ("Does radiation's effect differ between treatments?", "Control vs IFA", "Interaction: IFA vs Control", "", None),
    ("Does radiation's effect differ between treatments?", "Control vs Inosine", "Interaction: Inosine vs Control", "", None),
    ("Does radiation's effect differ between treatments?", "Inosine vs IFA", "Interaction: IFA vs Inosine", "", None),
]


def table():
    d = pd.read_csv(TABLE, sep="\t", keep_default_na=False, na_values=["NA"])
    for c in ("call", "robust"):
        d[c] = d[c].astype(str).str.upper().eq("TRUE")
    d["stratum"] = d.stratum.fillna("")
    return d


def rows():
    d = table()
    ab = d[(d.model == "abundance") & d.p.notna()]
    out = []
    for question, lab, row, st, swatch in ROWS:
        x = ab[(ab.row == row) & (ab.stratum == st)]
        h = x[x.call]
        up, down = h[h.coef > 0], h[h.coef < 0]
        out.append(dict(
            question=question, label=lab, row=row, stratum=st, swatch=swatch, tested=len(x),
            n=len(h), n_robust=int(h.robust.sum()),
            max_up=up.coef.max() if len(up) else None,
            max_down=-down.coef.min() if len(down) else None,
            min_q=h.q.min() if len(h) else None, typical=x.coef.median(),
            nominal=int((x.p < ALPHA).sum()),
            hits=h.sort_values("q")[["taxon", "coef", "q", "loo_holds_in", "robust"]].to_dict("records"),
            all=x))
    return out


def signed(v):
    return "0.00" if abs(v) < 0.005 else f"{v:+.2f}".replace("-", "−")


def html_rows(rs):
    """The table's <tr> lines, as the report shows them."""
    lines = []
    for i, r in enumerate(rs):
        sw = f'<span class="swatch" style="background:{r["swatch"]}"></span>' if r["swatch"] else ""
        ud = "—" if not r["n"] else (f'{signed(r["max_up"]) if r["max_up"] is not None else "—"} / '
                                    f'{signed(-r["max_down"]) if r["max_down"] is not None else "—"}')
        mq = "—" if r["min_q"] is None else f'{r["min_q"]:.1e}'
        chip = ('<span class="chip ok">robust</span>' if r["n_robust"] else
                '<span class="chip warn">not robust</span>' if r["n"] else
                '<span class="chip warn">below threshold</span>')
        cells = (f'<td>{sw}{r["label"]}</td><td class="n">{r["tested"]}</td><td class="n">{r["n"]} ({r["n_robust"]})</td>'
                 f'<td class="n">{ud}</td><td class="n">{mq}</td><td class="n">{signed(r["typical"])}</td><td>{chip}</td>')
        if i % 3 == 0:
            lines.append(f'<tr class="grp"><td rowspan="3" class="term q">{r["question"]}</td>{cells}</tr>')
        else:
            lines.append(f"<tr>{cells}</tr>")
    return lines


if __name__ == "__main__":
    rs = rows()
    for r in rs:
        print(f"{r['label']:22} {r['stratum']:13} tested={r['tested']} calls={r['n']} robust={r['n_robust']} "
              f"up={r['max_up']} down={r['max_down']} minq={r['min_q']} typical={r['typical']:.2f} nominal={r['nominal']}")
        for h in r["hits"]:
            print(f"      {h['taxon'].split(';')[-2:]} {h['coef']:+.2f} q={h['q']:.3g} loo={h['loo_holds_in']}/{N_MICE}")
    print("\n".join(html_rows(rs)))
