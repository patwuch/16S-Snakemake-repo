#!/usr/bin/env python3
"""Figures for the Day17/Day25 report's differential-abundance section, from MaAsLin3.

  1. Clustered heatmap: the 49 taxa seen in every Arm x Timepoint cell. Left, each cell's mean
     log2 share of reads over the samples that carry the taxon (what MaAsLin3's abundance model
     fits), centred per taxon; right, the seven report comparisons (abundance coefficient,
     log2), q < 0.05 filled dot, raw p < 0.05 ring. Rows clustered on the cell means.
  2. Forest plots, one per question, comparisons as facets, the heatmap's row order: coefficient
     with a 95% interval (coef +/- 1.96 SE, as MaAsLin3 reports it), and a dashed line at the
     median the test compares each coefficient against.

Each is drawn twice: all 49 taxa, and "_gplus" with only taxa named no further than genus (the
same fits, filtered, not refitted).

Usage, inside the qiime2 container:
    python3 plot_maaslin3_day17_day25.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage  # noqa: E402

from structural_zeros_day17_day25 import label  # noqa: E402

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-day17-day25")
TAB, FIG = os.path.join(OUT, "tables"), os.path.join(OUT, "figures")
STEM = "MaAsLin3_day17day25_Greengenes2"
ALPHA = 0.05

TEXT_INK, MUTED_INK, NULL_INK, RULE_INK = "#2b2b28", "#6b6b66", "#9a9a97", "#c9c9c4"
CELLS = ["IFA Day17", "IFA Day25", "IFA+IR Day17", "IFA+IR Day25"]
CELL_COLOR = {"IFA": "#90BFF9", "IFA+IR": "#0000FF"}
# Mark colours: darker steps of the arm hues, readable at marker size on white.
HUE = {"IFA": "#2f78d0", "IFA+IR": "#0a1fb0", "both": "#4d4d4a", "interaction": "#4a3aa7"}
DIVERGING = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f19c9b", "#b52d2c"])

# (key, question, [(facet label, report row, hue)], x label)
QUESTIONS = [
    ("time", "Does the community change from Day 17 to Day 25?",
     [("In IFA", "Day 17 to 25 in IFA", "IFA"),
      ("In IFA+IR", "Day 17 to 25 in IFA+IR", "IFA+IR"),
      ("Averaged over both arms", "Average Day 17 to 25 shift", "both")],
     "log2 coefficient, Day 25 vs Day 17"),
    ("arm", "Do the arms differ?",
     [("At Day 17", "IFA+IR vs IFA at Day 17", "IFA+IR"),
      ("At Day 25", "IFA+IR vs IFA at Day 25", "IFA+IR"),
      ("Averaged over both days", "Average IFA+IR vs IFA", "both")],
     "log2 coefficient, IFA+IR vs IFA"),
    ("interaction", "Does the Day 17 to 25 change differ between arms?",
     [("IFA+IR change minus IFA change", "Interaction (IFA+IR shift minus IFA shift)", "interaction")],
     "difference in the Day 25 vs Day 17 coefficient, IFA+IR minus IFA"),
]
SHORT = {"time": "Day 25 vs Day 17", "arm": "IFA+IR vs IFA", "interaction": "Interaction"}
# Heatmap column names, short enough to sit under the title.
HEAT_LABEL = {"Day 17 to 25 in IFA": "in IFA", "Day 17 to 25 in IFA+IR": "in IFA+IR",
              "Average Day 17 to 25 shift": "average", "IFA+IR vs IFA at Day 17": "at Day 17",
              "IFA+IR vs IFA at Day 25": "at Day 25", "Average IFA+IR vs IFA": "average",
              "Interaction (IFA+IR shift minus IFA shift)": "IFA+IR − IFA"}


def genus_or_above(taxon):
    r = [p.strip() for p in taxon.split(";") if p.strip() and not p.strip().endswith("__")]
    return not r[-1].startswith("s__")


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


res = pd.read_csv(os.path.join(TAB, f"{STEM}.tsv"), sep="\t")
res = res[res.model == "abundance"]
S1 = sorted(set(pd.read_csv(os.path.join(TAB, f"{STEM}_all_sets.tsv"), sep="\t").query("set == 'S1'").feature))
COLS = [(key, lab, row, hue) for key, _, comps, _ in QUESTIONS for lab, row, hue in comps]

# Cell means of log2 share where present, from the same counts MaAsLin3 fitted.
cnt = pd.read_csv(os.path.join(TAB, "collapsed_counts_day17day25_Greengenes2.tsv"), sep="\t", index_col=0)
smp = pd.read_csv(os.path.join(TAB, "samples_day17day25.tsv"), sep="\t", index_col=0).loc[cnt.columns]
cell = smp.Arm.map({"IFA": "IFA", "IFA_IR": "IFA+IR"}) + " " + smp.Timepoint
log_share = np.log2((cnt / cnt.sum()).where(cnt > 0))
cell_mean = pd.DataFrame({c: log_share.loc[:, cell == c].mean(axis=1) for c in CELLS}).loc[S1]


def coef_table(taxa, col):
    out = {}
    for key, lab, row, hue in COLS:
        d = res[res.row == row].set_index("feature").reindex(taxa)
        out[(key, lab)] = d[col].values
    return pd.DataFrame(out, index=taxa)


def heatmap(taxa, name, title):
    centred = cell_mean.loc[taxa].sub(cell_mean.loc[taxa].mean(axis=1), axis=0)
    Z = linkage(centred.values, method="average", metric="euclidean")
    order = [taxa[i] for i in leaves_list(Z)]
    centred = centred.loc[order]
    est, p, q = (coef_table(order, c) for c in ("coef", "pval_individual", "q"))
    skipped = coef_table(order, "skipped_median").fillna(False).astype(bool)
    n, lim = len(order), 4.0
    fig = plt.figure(figsize=(12, 2.8 + 0.23 * n))
    h = fig.get_figheight()
    top, bot = 1 - 1.75 / h, 1.15 / h
    ax_d = fig.add_axes([0.01, bot, 0.07, top - bot])
    ax_g = fig.add_axes([0.085, bot, 0.17, top - bot])
    ax_e = fig.add_axes([0.27, bot, 0.33, top - bot])
    dendrogram(Z, ax=ax_d, orientation="left", no_labels=True, color_threshold=0, above_threshold_color=MUTED_INK)
    ax_d.invert_yaxis()
    ax_d.axis("off")
    ax_g.imshow(centred.values.clip(-lim, lim), cmap=DIVERGING, vmin=-lim, vmax=lim, aspect="auto",
                interpolation="nearest")
    ax_g.set_xticks(range(len(CELLS)), CELLS, rotation=45, ha="left", fontsize=8.5)
    ax_g.xaxis.tick_top()
    ax_g.tick_params(axis="x", pad=16)
    for i, c in enumerate(CELLS):
        ax_g.add_patch(plt.Rectangle((i - 0.5, -1.15), 1, 0.55, color=CELL_COLOR[c.split()[0]], clip_on=False))
    ax_g.axvline(1.5, color="white", lw=2.5)
    ax_g.set_yticks([])
    ax_e.imshow(est.values.clip(-lim, lim), cmap=DIVERGING, vmin=-lim, vmax=lim, aspect="auto",
                interpolation="nearest")
    for (i, j), v in np.ndenumerate(q.values):
        if skipped.values[i, j]:
            continue
        if v < ALPHA:
            ax_e.plot(j, i, "o", ms=6.5, color=TEXT_INK, mec="white", mew=1.2)
        elif p.values[i, j] < ALPHA:
            ax_e.plot(j, i, "o", ms=5.5, mfc="none", mec=TEXT_INK, mew=1.1)
    ax_e.set_xticks(range(len(COLS)), [HEAT_LABEL[c[2]] for c in COLS], rotation=45, ha="left", fontsize=8.5)
    ax_e.xaxis.tick_top()
    for x in (2.5, 5.5):
        ax_e.axvline(x, color="white", lw=4)
    for k, (key, _, comps, _) in enumerate(QUESTIONS):
        ax_e.text(3 * k + (len(comps) - 1) / 2, n - 0.3, SHORT[key], ha="center", va="top", fontsize=8.5,
                  color=MUTED_INK)
    ax_e.set_yticks(range(n), [label(t) for t in order], fontsize=8.5)
    ax_e.yaxis.tick_right()
    for tl, t in zip(ax_e.get_yticklabels(), order):
        tl.set_fontstyle("italic")
        if (q.loc[t] < ALPHA).any():
            tl.set_fontweight("bold")
    for a in (ax_g, ax_e):
        a.tick_params(length=0, colors=TEXT_INK)
        for s in a.spines.values():
            s.set_visible(False)
    cax = fig.add_axes([0.085, bot - 0.75 / h, 0.17, 0.12 / h])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=DIVERGING, norm=plt.Normalize(-lim, lim)), cax=cax,
                      orientation="horizontal")
    cb.set_label(f"log2 (clipped at ±{lim:g})", fontsize=8.5, color=MUTED_INK)
    cb.ax.tick_params(labelsize=8, colors=MUTED_INK)
    cb.outline.set_visible(False)
    fig.text(0.27, bot - 0.72 / h, "●  q < 0.05 (BH within the comparison)      ○  raw p < 0.05, not a hit",
             fontsize=8.5, color=TEXT_INK, va="center")
    fig.suptitle(title, x=0.01, y=1 - 0.1 / h, ha="left", fontsize=11.5, color=TEXT_INK)
    fig.text(0.01, 1 - 0.42 / h, "Left: each group's mean log2 share of reads where present, centred per taxon.   "
             "Right: MaAsLin3 abundance coefficients for the seven comparisons.", ha="left", fontsize=9,
             color=MUTED_INK)
    save(fig, name)
    return order


def forest(order, key, question, comps, xlabel, name):
    n = len(order)
    fig, axes = plt.subplots(1, len(comps), figsize=(4.5 * len(comps) + (4.5 if len(comps) == 1 else 0),
                                                     1.6 + 0.23 * n), sharey=True, squeeze=False)
    axes = axes[0]
    lo_all, hi_all = [], []
    for ax, (lab, row, hue) in zip(axes, comps):
        d = res[res.row == row].set_index("feature").reindex(order)
        d["lo"], d["hi"] = d.coef - 1.96 * d.stderr, d.coef + 1.96 * d.stderr
        med = d.null_hypothesis.dropna().unique()
        ax.axvline(0, color=RULE_INK, lw=1, zorder=0)
        if len(med):
            ax.axvline(med[0], color=MUTED_INK, lw=1, ls=(0, (4, 3)), zorder=0)
        for i, (_, r) in enumerate(d.iterrows()):
            if pd.isna(r.coef):
                continue
            ok = not bool(r.skipped_median)
            hit, nom = ok and r.q < ALPHA, ok and r.pval_individual < ALPHA
            c = HUE[hue] if (hit or nom) else NULL_INK
            ax.plot([r.lo, r.hi], [i, i], color=c, lw=1.6 if (hit or nom) else 1, solid_capstyle="round")
            if hit:
                ax.plot(r.coef, i, "o", ms=7, color=c, mec="white", mew=1.2, zorder=3)
            elif nom:
                ax.plot(r.coef, i, "o", ms=6.5, mfc="white", mec=c, mew=1.6, zorder=3)
            else:
                ax.plot(r.coef, i, "o" if ok else "x", ms=4, color=c, zorder=3)
        lo_all += list(d.lo.dropna())
        hi_all += list(d.hi.dropna())
        ax.set_title(lab, fontsize=10, color=HUE[hue], loc="left", fontweight="bold")
        ax.set_yticks(range(n), [label(t) for t in order], fontsize=8.5)
        ax.grid(axis="x", color=RULE_INK, lw=0.5, alpha=0.6)
        ax.tick_params(length=0, colors=TEXT_INK, labelsize=8.5)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(RULE_INK)
    lim = np.ceil(max(abs(np.percentile(lo_all, 1)), abs(np.percentile(hi_all, 99))))
    for ax, (lab, row, hue) in zip(axes, comps):
        ax.set_xlim(-lim, lim)
        d = res[res.row == row].set_index("feature").reindex(order)
        for i, (_, r) in enumerate(d.iterrows()):
            if pd.isna(r.coef):
                continue
            if r.coef - 1.96 * r.stderr < -lim:
                ax.annotate("", (-lim, i), (-lim + 0.6, i), arrowprops=dict(arrowstyle="->", color=NULL_INK, lw=1))
            if r.coef + 1.96 * r.stderr > lim:
                ax.annotate("", (lim, i), (lim - 0.6, i), arrowprops=dict(arrowstyle="->", color=NULL_INK, lw=1))
    axes[0].invert_yaxis()
    for tl, t in zip(axes[0].get_yticklabels(), order):
        tl.set_fontstyle("italic")
        hits = res[res.row.isin([c[1] for c in comps]) & (res.feature == t) & ~res.skipped_median.astype(bool)]
        if (hits.q < ALPHA).any():
            tl.set_fontweight("bold")
    fig.supxlabel(xlabel, fontsize=9.5, color=MUTED_INK)
    h = fig.get_figheight()
    fig.suptitle(question, x=0.01, y=1 - 0.05 / h, ha="left", va="top", fontsize=11.5, color=TEXT_INK)
    fig.text(0.01, 1 - 0.38 / h, "●  q < 0.05 (BH within the comparison)   ○  raw p < 0.05   ×  not tested against the median\n"
             "bars: 95% interval   dashed line: the median each coefficient is tested against",
             ha="left", va="top", fontsize=8.5, color=MUTED_INK)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.85 / h))
    save(fig, name)


if __name__ == "__main__":
    for suffix, taxa, what in (("", S1, f"the {len(S1)} taxa seen in every arm × day group"),
                               ("_gplus", [t for t in S1 if genus_or_above(t)], "taxa named no further than genus")):
        order = heatmap(taxa, f"{STEM}{suffix}_heatmap", f"MaAsLin3 abundance by group and by comparison: {what}")
        for key, question, comps, xlabel in QUESTIONS:
            forest(order, key, question, comps, xlabel, f"{STEM}{suffix}_forest_{key}")
