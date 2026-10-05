#!/usr/bin/env python3
"""Targeted PICRUSt2 panels for ahmad260430: inosine and isoferulic acid (picrust2_targets.py).

For each report and compound, one grid: the curated features (rows, grouped by database and role)
against the report's MaAsLin3 rows (columns), from run_picrust2_maaslin3.R's tables. Colour is the
change in share of the predicted total (log2, abundance model); a filled ring marks a robust call
(q < 0.05 in every leave-one-mouse-out refit), a hollow ring a call that is not robust, and a grey
cell a feature not tested in that row. Each database block opens with the row's typical shift, the
median change over all its tested features, so a feature can be read against the whole profile.

Writes <report>/figures/picrust2_<compound>.{png,svg} and <report>/tables/PICRUSt2_<compound>_results.tsv.

Usage, inside the qiime2 container:
    python3 plot_picrust2_targets.py factorial|day17-day25
"""
import html
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from picrust2_targets import FEATURES, TITLES  # noqa: E402

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(THIS_DIR)
REPORTS = {
    "factorial": ("comparison-Treatment-x-Radiation-factorial", [
        ("IFA vs Control", "None", "IFA vs Control"), ("Inosine vs Control", "None", "Inosine vs Control"),
        ("IFA vs Inosine", "None", "IFA vs Inosine"), ("IFA vs Control", "IR", "IFA+IR vs IR"),
        ("Inosine vs Control", "IR", "Inosine+IR vs IR"), ("IFA vs Inosine", "IR", "IFA+IR vs Inosine+IR"),
        ("IR within Control", "", "IR vs Control"), ("IR within IFA", "", "IFA+IR vs IFA"),
        ("IR within Inosine", "", "Inosine+IR vs Inosine"), ("Interaction: IFA vs Control", "", "Interaction: IFA vs Control"),
        ("Interaction: Inosine vs Control", "", "Interaction: Inosine vs Control"),
        ("Interaction: IFA vs Inosine", "", "Interaction: IFA vs Inosine")]),
    "day17-day25": ("comparison-day17-day25", [
        ("Day 17 to 25 in IFA", "", "Day 17 to 25, IFA"), ("Day 17 to 25 in IFA+IR", "", "Day 17 to 25, IFA+IR"),
        ("IFA+IR vs IFA at Day 17", "", "IFA+IR vs IFA, Day 17"), ("IFA+IR vs IFA at Day 25", "", "IFA+IR vs IFA, Day 25"),
        ("Interaction (IFA+IR shift minus IFA shift)", "", "Interaction (shift difference)"),
        ("Average Day 17 to 25 shift", "", "Average Day 17 to 25 shift"), ("Average IFA+IR vs IFA", "", "Average IFA+IR vs IFA")]),
}
DBS = ["MetaCyc", "EC", "KO"]
DBNAME = {"MetaCyc": "MetaCyc pathways", "EC": "Enzymes (EC)", "KO": "KEGG orthologs (KO)"}
INK, MUTED, NA = "#14181d", "#6b7280", "#e5e5e2"


def load(report_dir):
    out = {}
    for lv in DBS:
        d = pd.read_csv(os.path.join(STUDY, report_dir, "tables", f"PICRUSt2_{lv}_maaslin3_rows.tsv"), sep="\t",
                        keep_default_na=False, na_values=["NA"])
        out[lv] = d[d.model == "abundance"].copy()
    return out


def main(report):
    rdir, cols = REPORTS[report]
    res = load(rdir)
    med = {(lv, r, s): res[lv][(res[lv].row == r) & (res[lv].stratum == s) & res[lv].p.notna()].coef.median()
           for lv in DBS for r, s, _ in cols}
    for comp, feats in FEATURES.items():
        rows, recs = [], []   # rows: (kind, db, label) with kind in {"db", "median", "role", "feat"}
        for lv in DBS:
            fl = [f for f in feats if f[0] == lv]
            if not fl:
                continue
            rows.append(("median", lv, f"{DBNAME[lv]}: typical shift"))
            role = None
            for _db, fid, rl, lab in fl:
                if rl != role:
                    rows.append(("role", lv, rl)); role = rl
                rows.append(("feat", lv, (fid, lab)))
                for r, s, clab in cols:
                    x = res[lv][(res[lv].row == r) & (res[lv].stratum == s) & (res[lv].id == fid)]
                    x = x.iloc[0] if len(x) else None
                    tested = x is not None and pd.notna(x.p)
                    recs.append(dict(compound=comp, database=lv, feature=fid, role=rl, label=lab,
                                     description=html.unescape(str(x.description)) if x is not None else "",
                                     row=clab, coef=x.coef if tested else np.nan, p=x.p if tested else np.nan,
                                     q=x.q if tested else np.nan, call=bool(tested and x.call), robust=bool(tested and x.robust),
                                     loo_holds_in=x.loo_holds_in if tested else np.nan,
                                     n_predicted=x.n_predicted if x is not None else 0, typical_shift=med[(lv, r, s)]))
        tab = pd.DataFrame(recs)
        tab.to_csv(os.path.join(STUDY, rdir, "tables", f"PICRUSt2_{comp}_results.tsv"), sep="\t", index=False)

        vals = tab.coef.dropna().abs()
        lim = max(1.0, float(np.ceil(np.nanpercentile(vals, 97) if len(vals) else 1)))
        norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
        cmap = plt.get_cmap("RdBu_r")
        nr, nc = len(rows), len(cols)
        fig_h = 1.6 + 0.24 * nr
        fig = plt.figure(figsize=(3.6 + 0.55 * nc + 3.2, fig_h))
        ax = fig.add_axes([0.40, 0.9 / fig_h, 0.52 * nc / (nc + 4.5) + 0.02, 1 - 2.4 / fig_h])
        for i, (kind, lv, lab) in enumerate(rows):
            y = nr - 1 - i
            if kind == "role":
                ax.text(-0.6, y, lab, ha="right", va="center", fontsize=7.5, color=MUTED, style="italic")
                continue
            for j, (r, s, clab) in enumerate(cols):
                if kind == "median":
                    v = med[(lv, r, s)]
                    ax.add_patch(plt.Rectangle((j - 0.5, y - 0.5), 1, 1, color=cmap(norm(np.clip(v, -lim, lim))), lw=0))
                    ax.text(j, y, f"{v:+.1f}", ha="center", va="center", fontsize=6.5,
                            color="white" if abs(v) > 0.6 * lim else INK)
                    continue
                t = tab[(tab.database == lv) & (tab.feature == lab[0]) & (tab.row == clab)].iloc[0]
                if np.isnan(t.coef):
                    ax.add_patch(plt.Rectangle((j - 0.5, y - 0.5), 1, 1, color=NA, lw=0))
                    continue
                c = cmap(norm(np.clip(t.coef, -lim, lim)))
                ax.add_patch(plt.Rectangle((j - 0.5, y - 0.5), 1, 1, color=c, lw=0))
                # Marker colour contrasts with the cell: white on dark cells, ink on light ones.
                mk = "white" if 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2] < 0.45 else INK
                if t.robust:
                    ax.scatter(j, y, s=30, facecolor=mk, edgecolor="none", zorder=3)
                elif t.call:
                    ax.scatter(j, y, s=30, facecolor="none", edgecolor=mk, linewidth=1.1, zorder=3)
            if kind == "median":
                ax.axhline(y + 0.5, color=INK, lw=0.8)
                ax.text(-0.6, y, lab, ha="right", va="center", fontsize=8, color=INK, weight="bold")
            else:
                ax.text(-0.6, y, lab[1], ha="right", va="center", fontsize=7.5, color=INK)
        ax.set_xlim(-0.5, nc - 0.5)
        ax.set_ylim(-0.5, nr - 0.5)
        ax.set_xticks(range(nc), [c[2] for c in cols], rotation=55, ha="left", fontsize=7.5)
        ax.xaxis.tick_top()
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)
        cax = fig.add_axes([0.40, 0.35 / fig_h, 0.22, 0.12 / fig_h])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
        cb.set_label("change in share of the predicted total (log2; first-named group higher → red)", fontsize=7.5, color=MUTED)
        cb.ax.tick_params(labelsize=7)
        fig.text(0.66, 0.42 / fig_h, "●  robust call    ○  call, not robust    grey: not tested", fontsize=7.5, color=MUTED, va="center")
        fig.suptitle(f"{TITLES[comp]}: curated predicted features", x=0.01, ha="left", fontsize=11, color=INK)
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(STUDY, rdir, "figures", f"picrust2_{comp}.{ext}"), dpi=200, bbox_inches="tight",
                        facecolor="white")
        plt.close(fig)
        print(f"{report} {comp}: {len(tab.feature.unique())} features; robust calls {int(tab.robust.sum())}, "
              f"calls {int(tab.call.sum())}, tested cells {int(tab.coef.notna().sum())} of {len(tab)}")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in REPORTS:
        sys.exit("usage: plot_picrust2_targets.py factorial|day17-day25")
    main(sys.argv[1])
