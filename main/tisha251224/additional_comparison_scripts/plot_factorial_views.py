#!/usr/bin/env python3
"""Whole-design views of the MaAsLin3 abundance results: a clustered factorial
heatmap and a faceted forest plot, over every taxon with a call in any row.

Heatmap: one row per taxon, one column per Diet x Condition group, laid out as
the design (three diets, unirradiated then irradiated). A cell is the mean log2
share of reads over the group's samples where the taxon is present (what the
abundance model fits), centred on the taxon's own average, so a row reads as
where the taxon is relatively high or low. A group with the taxon in no sample
is hatched gray with "0/n"; one with it in only some samples prints "k/n".
Rows are clustered (average linkage, Euclidean, optimal leaf order); the strip
on the right marks, for each of the 12 tests, a robust call (filled) or a call
that does not hold when some one mouse is left out (hollow).

Forest: the same taxa in the same order, one facet per question, the three
contrasts of each question side by side in a row, each estimate with its 95%
interval (t, with each taxon's df recovered from its t statistic and p, as in
plot_interaction_forest.py). Filled = robust call, hollow = call, faint = no
call; a missing mark means the taxon was not testable in that row.

Usage, inside the qiime2 container:
    python3 plot_factorial_views.py [--genus-names-only]
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import t as t_dist

from da_rows import species_named
from plot_interaction_forest import implied_df

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
MA = os.path.join(OUT, "tables", "maaslin3")
FIG = os.path.join(OUT, "figures")
PRE = "MaAsLin3_factorial_Greengenes2"

DIETS, CONDS = ["Cellulose", "Inulin", "RS"], ["Control", "IR"]
INK, MUTED, RULE = "#2b2b28", "#6b6b66", "#c9c9c4"
ABSENT = "#d9d9d6"
# Diverging pair of the reference palette: blue (lower) <-> gray <-> red (higher).
CMAP = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])
# Contrast hues: the dark (irradiated) step of the first-named diet, a violet for
# RS vs Inulin, and each diet's own hue for irradiation within it. Same as the
# interaction forest plot.
HUE = {"Inulin": "#1f5bb0", "RS": "#b8460f", "RSvI": "#4a3aa7", "Cellulose": "#12805a"}
# (facet title, [(run_maaslin3_rows row, stratum, short label, hue)])
FACETS = [
    ("Diets, unirradiated", [("Inulin vs Cellulose", "Control", "Inulin vs Cellulose", HUE["Inulin"]),
                             ("RS vs Cellulose", "Control", "RS vs Cellulose", HUE["RS"]),
                             ("RS vs Inulin", "Control", "RS vs Inulin", HUE["RSvI"])]),
    ("Diets, irradiated", [("Inulin vs Cellulose", "IR", "Inulin+IR vs Cellulose+IR", HUE["Inulin"]),
                           ("RS vs Cellulose", "IR", "RS+IR vs Cellulose+IR", HUE["RS"]),
                           ("RS vs Inulin", "IR", "RS+IR vs Inulin+IR", HUE["RSvI"])]),
    ("Irradiation within diet", [("IR within Cellulose", "", "Cellulose+IR vs Cellulose", HUE["Cellulose"]),
                                 ("IR within Inulin", "", "Inulin+IR vs Inulin", HUE["Inulin"]),
                                 ("IR within RS", "", "RS+IR vs RS", HUE["RS"])]),
    ("Interaction: diet vs diet",
     [("Interaction: Inulin vs Cellulose", "", "Inulin vs Cellulose", HUE["Inulin"]),
      ("Interaction: RS vs Cellulose", "", "RS vs Cellulose", HUE["RS"]),
      ("Interaction: RS vs Inulin", "", "RS vs Inulin", HUE["RSvI"])]),
]
STRIP = ["I–C", "R–C", "R–I", "I–C", "R–C", "R–I", "C", "I", "R", "I–C", "R–C", "R–I"]


def label(taxon):
    """Deepest assigned rank, named by that rank ("genus X (no species call)")."""
    parts = [p for p in taxon.split(";") if not p.endswith("__")]
    rank, name = parts[-1][0], re.sub(r"_\d+(?=$| )", "", parts[-1][3:])
    if rank == "s":
        return name.replace("_", " ")
    below = {"g": "species", "f": "genus", "o": "family"}[rank]
    word = {"g": "genus", "f": "family", "o": "order"}[rank]
    return f"{word} {name.replace('_', ' ')} (no {below} call)"


def load():
    res = pd.read_csv(os.path.join(MA, f"{PRE}_min2_rows.tsv"), sep="\t")
    res["stratum"] = res.stratum.fillna("")
    res = res[res.model == "abundance"]
    counts = pd.read_csv(os.path.join(MA, "counts_Greengenes2.tsv"), sep="\t", index_col=0)
    counts = counts.loc[:, (counts > 0).sum() >= 2]          # the features MaAsLin3 was given
    meta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t", index_col=0,
                       keep_default_na=False).loc[counts.index]
    return res, counts, meta


def cell_matrix(counts, meta, feats):
    share = counts.div(counts.sum(axis=1), axis=0)
    cols, mean, pres = [], {}, {}
    for c in CONDS:
        for d in DIETS:
            s = meta.index[(meta.Diet == d) & (meta.Condition == c)]
            g = meta.loc[s[0], "Group"]
            cols.append(g)
            x = share.loc[s, feats]
            mean[g] = np.log2(x.where(x > 0)).mean()
            pres[g] = [(int((x[f] > 0).sum()), len(s)) for f in feats]
    m = pd.DataFrame(mean)[cols]
    return m.sub(m.mean(axis=1), axis=0), pd.DataFrame(pres, index=feats)[cols], cols


def draw_heatmap(M, P, cols, marks, labels, colors, path, title):
    n = len(M)
    fig = plt.figure(figsize=(12.5, 0.30 * n + 2.6))
    H = fig.get_figheight()
    top, bot = 1 - 1.45 / H, 0.85 / H
    ax_d = fig.add_axes([0.02, bot, 0.10, top - bot])
    ax_h = fig.add_axes([0.125, bot, 0.30, top - bot])
    ax_s = fig.add_axes([0.445, bot, 0.21, top - bot])
    ax_c = fig.add_axes([0.125, 0.30 / H, 0.30, 0.13 / H])

    Z = linkage(M.fillna(0).values, method="average", metric="euclidean", optimal_ordering=True)
    dn = dendrogram(Z, orientation="left", ax=ax_d, no_labels=True, color_threshold=0,
                    above_threshold_color=MUTED, link_color_func=lambda k: MUTED)
    order = dn["leaves"][::-1]                                 # top row first
    ax_d.set_ylim(0, 10 * n); ax_d.invert_yaxis(); ax_d.axis("off")

    lim = np.nanpercentile(np.abs(M.values), 97)
    norm = TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim)
    xpos = [0, 1, 2, 3.25, 4.25, 5.25]                         # gap between the two conditions
    for i, r in enumerate(order):
        for j, g in enumerate(cols):
            v, (k, tot) = M.iloc[r][g], P.iloc[r][g]
            if k == 0:
                ax_h.add_patch(Rectangle((xpos[j], i), 1, 1, facecolor=ABSENT, edgecolor="white", lw=1.2, hatch="///"))
                ax_h.text(xpos[j] + 0.5, i + 0.5, f"0/{tot}", ha="center", va="center", fontsize=6.5, color=INK)
            else:
                ax_h.add_patch(Rectangle((xpos[j], i), 1, 1, facecolor=CMAP(norm(v)), edgecolor="white", lw=1.2))
                if k < tot:
                    ax_h.text(xpos[j] + 0.5, i + 0.5, f"{k}/{tot}", ha="center", va="center", fontsize=6.5, color=INK)
    ax_h.set_xlim(0, 6.25); ax_h.set_ylim(n, 0); ax_h.axis("off")
    for j, g in enumerate(cols):
        ax_h.add_patch(Rectangle((xpos[j] + 0.1, -0.55), 0.8, 0.28, facecolor=colors[g], edgecolor="none", clip_on=False))
        ax_h.text(xpos[j] + 0.5, -0.75, g.replace("+IR", ""), ha="center", va="bottom", fontsize=8, color=INK, rotation=0)
    ax_h.text(1.5, -1.55, "Unirradiated", ha="center", fontsize=9, color=INK, weight="bold")
    ax_h.text(4.75, -1.55, "Irradiated", ha="center", fontsize=9, color=INK, weight="bold")

    # Test strip: 12 columns in four groups of three.
    sx = [0, 1, 2, 3.4, 4.4, 5.4, 6.8, 7.8, 8.8, 10.2, 11.2, 12.2]
    for i, r in enumerate(order):
        for j, x in enumerate(sx):
            state = marks[j].get(M.index[r])
            if state == "robust":
                ax_s.plot(x + 0.5, i + 0.5, "o", ms=6.5, color=INK)
            elif state == "call":
                ax_s.plot(x + 0.5, i + 0.5, "o", ms=6.5, mfc="white", mec=INK, mew=1.2)
            elif state is None:
                ax_s.plot(x + 0.5, i + 0.5, "_", ms=5, color=RULE)
        ax_s.text(13.6, i + 0.5, labels[M.index[r]], ha="left", va="center", fontsize=8, color=INK, style="italic")
    for gi, (t, _) in enumerate(FACETS):
        x0 = sx[3 * gi]
        ax_s.add_patch(Rectangle((x0 - 0.1, 0), 3.2, n, facecolor="#f4f4f2", edgecolor="none", zorder=-1))
        ax_s.text(x0 + 1.5, -1.55, ["Diets\nunirr.", "Diets\nirr.", "Irradiation\nwithin", "Inter-\naction"][gi],
                  ha="center", fontsize=7.5, color=INK, weight="bold")
    for j, x in enumerate(sx):
        ax_s.text(x + 0.5, -0.4, STRIP[j], ha="center", va="bottom", fontsize=6.5, color=MUTED, rotation=90)
    ax_s.set_xlim(-0.2, 13.4); ax_s.set_ylim(n, 0); ax_s.axis("off")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cb = fig.colorbar(sm, cax=ax_c, orientation="horizontal")
    cb.set_label("mean log2 share of reads where present, relative to the taxon's own average", fontsize=7.5, color=INK)
    cb.ax.tick_params(labelsize=7, colors=MUTED); cb.outline.set_edgecolor(RULE)
    fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=6.5, color=INK, label="robust call"),
                        Line2D([], [], ls="", marker="o", ms=6.5, mfc="white", mec=INK, label="call, lost when one mouse is left out"),
                        Line2D([], [], ls="", marker="_", ms=6, color=RULE, label="not testable (absent from a compared group)"),
                        Rectangle((0, 0), 1, 1, facecolor=ABSENT, hatch="///", edgecolor="white", label="absent from the group"),
                        Line2D([], [], ls="", label="blank: tested, no call")],
               loc="lower left", bbox_to_anchor=(0.445, 0.05 / H * 2), ncol=2, frameon=False, fontsize=7.5)
    fig.suptitle(title, x=0.02, ha="left", y=1 - 0.22 / H, fontsize=12, color=INK)
    fig.text(0.02, 1 - 0.50 / H, f"MaAsLin3 · Greengenes2 · {n} taxa with a call at q < 0.05 in any row · rows clustered",
             fontsize=8, color=MUTED)
    for ext in ("png", "svg"):
        fig.savefig(f"{path}.{ext}", dpi=200, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    return [M.index[r] for r in order]


def draw_forest(res, order, labels, path, title):
    n = len(order)
    fig, axes = plt.subplots(1, 4, sharey=True, figsize=(15, 0.36 * n + 2.4))
    H = fig.get_figheight()
    y = {f: i for i, f in enumerate(order)}
    for ax, (ftitle, cons) in zip(axes, FACETS):
        lo_hi = []
        for off, (row, st, short, hue) in zip((-0.25, 0.0, 0.25), cons):
            d = res[(res.row == row) & (res.stratum == st) & res.feature.isin(order) & res.p.notna()]
            for r in d.itertuples():
                df = implied_df(r.coef / r.stderr, r.p)
                half = t_dist.ppf(0.975, df) * r.stderr
                yy = y[r.feature] + off
                faint = not r.call
                ax.plot([r.coef - half, r.coef + half], [yy, yy], color=hue, lw=1.4, alpha=0.3 if faint else 1, solid_capstyle="round")
                if r.robust:
                    ax.plot(r.coef, yy, "o", ms=6, color=hue, mec="white", mew=0.8, zorder=3)
                elif r.call:
                    ax.plot(r.coef, yy, "o", ms=6, mfc="white", mec=hue, mew=1.5, zorder=3)
                else:
                    ax.plot(r.coef, yy, "o", ms=4, color=hue, alpha=0.35, mew=0, zorder=2)
                if df >= 3:
                    lo_hi += [r.coef - half, r.coef + half]
        m = max(abs(v) for v in lo_hi) if lo_hi else 1
        m = min(m, 14)
        ax.set_xlim(-m * 1.05, m * 1.05)
        ax.axvline(0, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=1)
        for i in range(n):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color="#f4f4f2", zorder=0, lw=0)
        ax.set_title(ftitle, fontsize=9.5, color=INK, loc="left", pad=50, weight="bold")
        ax.legend(handles=[Line2D([], [], ls="-", marker="o", color=h, ms=5, lw=1.4, label=s) for _, _, s, h in cons],
                  loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, frameon=False, fontsize=7, handlelength=1.6,
                  borderaxespad=0.2, labelcolor=INK)
        ax.tick_params(axis="x", labelsize=7.5, colors=MUTED)
        ax.grid(axis="x", color=RULE, lw=0.5)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color(RULE)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(range(n))
    axes[0].set_yticklabels([labels[f] for f in order], fontsize=8, color=INK, style="italic")
    axes[0].set_ylim(n - 0.5, -0.5)
    fig.supxlabel("log2 change in share of reads, first-named relative to second (95% CI)", fontsize=9, color=INK, y=0.25 / H)
    fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=6, color=MUTED, label="robust call"),
                        Line2D([], [], ls="", marker="o", ms=6, mfc="white", mec=MUTED, mew=1.5, label="call, lost when one mouse is left out"),
                        Line2D([], [], ls="", marker="o", ms=4, color=MUTED, alpha=0.35, label="no call")],
               loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=3, frameon=False, fontsize=8)
    fig.suptitle(title, x=0.01, ha="left", fontsize=12, color=INK, y=1 - 0.15 / H)
    fig.subplots_adjust(left=0.2, right=0.99, wspace=0.12, top=1 - 1.25 / H, bottom=1.0 / H)
    for ext in ("png", "svg"):
        fig.savefig(f"{path}.{ext}", dpi=200, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genus-names-only", action="store_true",
                    help="same fit and calls; keep only taxa named no further than genus")
    a = ap.parse_args()
    res, counts, meta = load()
    feats = sorted(res.loc[res.call, "feature"].unique())
    key = pd.read_csv(os.path.join(MA, "feature_key.tsv"), sep="\t").set_index("feature").taxon
    if a.genus_names_only:
        feats = [f for f in feats if not species_named(key[f])]
    labels = {f: label(key[f]) for f in feats}
    M, P, cols = cell_matrix(counts, meta, feats)
    colors = dict(zip(meta.Group, meta.PrimaryColor))
    marks = []
    for _, cons in FACETS:
        for row, st, _s, _h in cons:
            d = res[(res.row == row) & (res.stratum == st) & res.feature.isin(feats)].set_index("feature")
            marks.append({f: ("robust" if d.loc[f, "robust"] else "call" if d.loc[f, "call"] else
                              ("none" if pd.notna(d.loc[f, "p"]) else None)) if f in d.index else None for f in feats})
    tag = "_gplus" if a.genus_names_only else ""
    view = " · genus-level names only" if a.genus_names_only else ""
    order = draw_heatmap(M, P, cols, marks, labels, colors, os.path.join(FIG, f"{PRE}{tag}_factorial_heatmap"),
                         f"Where each taxon's share is high or low across the design{view}")
    draw_forest(res, order, labels, os.path.join(FIG, f"{PRE}{tag}_factorial_forest"),
                f"Every test of every taxon with a call{view}")
    print(f"  {len(feats)} taxa -> {PRE}{tag}_factorial_heatmap / _factorial_forest")


if __name__ == "__main__":
    main()
