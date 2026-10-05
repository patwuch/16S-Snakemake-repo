#!/usr/bin/env python3
"""Figures for the factorial report's differential-abundance section, from ALDEx3.

  1. Clustered factorial heatmap (abundance): the 46 tested taxa, left the
     ALDEx3 model-fitted mean of each group (log2, centred per taxon), right
     the twelve report comparisons (log2FC), q < 0.05 filled dot, raw p < 0.05
     ring. Rows are clustered on the group means.
  2. Faceted forest plots, one per question: the three comparisons that answer
     it as facets, every tested taxon in the heatmap's row order, 95% intervals
     (quantiles of the draw mixture, so an interval excludes zero exactly when
     raw p < 0.05), filled dot q < 0.05 within the question's family.
  3. Clustered presence heatmap: every taxon seen in >= 6 mice (tested and
     untestable), mice with the taxon out of four per group. Taxa with a
     complete-separation call are marked, untestable ones flagged.

Each is drawn twice: all taxa, and "_gplus" with only taxa named no further
than genus (same fits, a filter, not a refit).

Usage, inside the qiime2 container:
    python3 plot_aldex3_factorial.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, LogNorm  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage  # noqa: E402
from scipy.stats import t as t_dist  # noqa: E402

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(THIS_DIR)
OUT = os.path.join(STUDY, "comparison-Treatment-x-Radiation-factorial")
TAB, FIG = os.path.join(OUT, "tables"), os.path.join(OUT, "figures")
TAG = "Greengenes2_Day17_minpres6"
STEM = f"ALDEx3_factorial_{TAG}_gamma0.5"
ALPHA = 0.05

TEXT_INK, MUTED_INK, NULL_INK, RULE_INK = "#2b2b28", "#6b6b66", "#9a9a97", "#c9c9c4"
# Darker steps of the study's treatment hues that pass the palette validator at
# mark size (as plot_interaction_forest.py); Control has no hue of its own.
HUE = {"IFA": "#1F6FD0", "Inosine": "#C25200", "Control": "#4d4d4a", "IFA vs Inosine": "#4a3aa7"}
GROUP_COLOR = {"Control": "#808080", "IR": "#000000", "IFA": "#90BFF9", "IFA+IR": "#0000FF",
               "Inosine": "#FF8000", "Inosine+IR": "#FF0000"}
GROUPS = ["Control", "IR", "Inosine", "Inosine+IR", "IFA", "IFA+IR"]  # radiation pairs adjacent
DIVERGING = LinearSegmentedColormap.from_list(
    "div", ["#1c5cab", "#86b6ef", "#f0efec", "#f19c9b", "#b52d2c"])

# The report's questions, each a family of three comparisons (label, row, stratum, hue).
QUESTIONS = [
    ("treatment_unirradiated", "Do the treatments differ in unirradiated mice?",
     [("IFA vs Control", "IFA vs Control", "unirradiated", "IFA"),
      ("Inosine vs Control", "Inosine vs Control", "unirradiated", "Inosine"),
      ("IFA vs Inosine", "IFA vs Inosine", "unirradiated", "IFA vs Inosine")]),
    ("treatment_irradiated", "Do the treatments differ in irradiated mice?",
     [("IFA+IR vs IR", "IFA vs Control", "irradiated", "IFA"),
      ("Inosine+IR vs IR", "Inosine vs Control", "irradiated", "Inosine"),
      ("IFA+IR vs Inosine+IR", "IFA vs Inosine", "irradiated", "IFA vs Inosine")]),
    ("radiation", "Does radiation change the microbiome within each treatment?",
     [("Within Control", "Radiation within Control", "", "Control"),
      ("Within IFA", "Radiation within IFA", "", "IFA"),
      ("Within Inosine", "Radiation within Inosine", "", "Inosine")]),
    ("interaction", "Does radiation's effect differ between treatments?",
     [("IFA vs Control", "Interaction: IFA vs Control", "", "IFA"),
      ("Inosine vs Control", "Interaction: Inosine vs Control", "", "Inosine"),
      ("IFA vs Inosine", "Interaction: IFA vs Inosine", "", "IFA vs Inosine")]),
]
X_LABEL = {"treatment_unirradiated": "log2 fold change, first-named vs second",
           "treatment_irradiated": "log2 fold change, first-named vs second",
           "radiation": "log2 fold change, irradiated vs unirradiated",
           "interaction": "difference in radiation's log2FC, first-named minus second"}


def ranks(taxon):
    return [p.strip() for p in taxon.split(";") if p.strip() and not p.strip().endswith("__")]


def genus_or_above(taxon):
    return not ranks(taxon)[-1].startswith("s__")


def label(taxon):
    """Deepest rank actually assigned; empty-slot paths say so."""
    r = ranks(taxon)
    name = r[-1][3:].replace("_", " ") if r[-1].startswith("s__") else r[-1][3:]
    if r[-1].startswith("s__"):
        return name
    empty_slot = taxon.rstrip().endswith("s__")
    if r[-1].startswith("g__"):
        return f"{name} ({'species unresolved' if empty_slot else 'genus only'})"
    rank = {"f": "family", "o": "order"}.get(r[-1][0], "higher rank")
    return f"{name} ({rank} only{', species unresolved' if empty_slot else ''})"


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


res = pd.read_csv(os.path.join(TAB, f"{STEM}.tsv"), sep="\t", keep_default_na=False)
res["p"], res["q"] = res.p.astype(float), res.q.astype(float)
gm = pd.read_csv(os.path.join(TAB, f"{STEM}_groupmeans.tsv"), sep="\t").set_index("taxon")[GROUPS]
COLS = [(lab, row, st, hue, key) for key, _, comps in QUESTIONS for lab, row, st, hue in comps]


def effect_matrix(taxa, col):
    out = {}
    for lab, row, st, hue, key in COLS:
        d = res[(res.row == row) & (res.stratum == st)].set_index("taxon").loc[taxa]
        out[(key, lab)] = d[col].values
    return pd.DataFrame(out, index=taxa)


def row_order(taxa):
    centred = gm.loc[taxa].sub(gm.loc[taxa].mean(axis=1), axis=0)
    Z = linkage(centred.values, method="average", metric="euclidean")
    return Z, [taxa[i] for i in leaves_list(Z)]


# ---------------------------------------------------------------- heatmap
def abundance_heatmap(taxa, name, title):
    Z, order = row_order(taxa)
    centred = gm.loc[order].sub(gm.loc[order].mean(axis=1), axis=0)
    est, p, q = (effect_matrix(order, c) for c in ("log2FC", "p", "q"))
    n = len(order)
    lim = 4.0
    fig = plt.figure(figsize=(13.5, 2.6 + 0.23 * n))
    top = 1 - 1.35 / fig.get_figheight()
    bot = 1.15 / fig.get_figheight()
    ax_d = fig.add_axes([0.01, bot, 0.07, top - bot])
    ax_g = fig.add_axes([0.085, bot, 0.20, top - bot])
    ax_e = fig.add_axes([0.30, bot, 0.43, top - bot])
    dendrogram(Z, ax=ax_d, orientation="left", no_labels=True, color_threshold=0,
               above_threshold_color=MUTED_INK, labels=taxa)
    ax_d.invert_yaxis()
    ax_d.axis("off")
    ax_g.imshow(centred.values.clip(-lim, lim), cmap=DIVERGING, vmin=-lim, vmax=lim, aspect="auto",
                interpolation="nearest")
    ax_g.set_xticks(range(len(GROUPS)), GROUPS, rotation=45, ha="left", fontsize=8.5)
    ax_g.xaxis.tick_top()
    ax_g.tick_params(axis="x", pad=16)
    for i, g in enumerate(GROUPS):
        ax_g.add_patch(plt.Rectangle((i - 0.5, -1.15), 1, 0.55, color=GROUP_COLOR[g], clip_on=False))
    for x in (1.5, 3.5):
        ax_g.axvline(x, color="white", lw=2.5)
    ax_g.set_yticks([])
    # Effects: 12 columns in four blocks of three.
    ax_e.imshow(est.values.clip(-lim, lim), cmap=DIVERGING, vmin=-lim, vmax=lim, aspect="auto",
                interpolation="nearest")
    for (i, j), v in np.ndenumerate(q.values):
        if v < ALPHA:
            ax_e.plot(j, i, "o", ms=6.5, color=TEXT_INK, mec="white", mew=1.2)
        elif p.values[i, j] < ALPHA:
            ax_e.plot(j, i, "o", ms=5.5, mfc="none", mec=TEXT_INK, mew=1.1)
    ax_e.set_xticks(range(len(COLS)), [c[0] for c in COLS], rotation=45, ha="left", fontsize=8.5)
    ax_e.xaxis.tick_top()
    for x in (2.5, 5.5, 8.5):
        ax_e.axvline(x, color="white", lw=4)
    for k, (key, qtext, _) in enumerate(QUESTIONS):
        short = {"treatment_unirradiated": "Treatment, unirradiated", "treatment_irradiated": "Treatment, irradiated",
                 "radiation": "Radiation, within", "interaction": "Interaction"}[key]
        ax_e.text(3 * k + 1, n - 0.3, short, ha="center", va="top", fontsize=8.5, color=MUTED_INK)
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
    cax = fig.add_axes([0.085, bot - 0.75 / fig.get_figheight(), 0.20, 0.12 / fig.get_figheight()])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=DIVERGING, norm=plt.Normalize(-lim, lim)), cax=cax,
                      orientation="horizontal")
    cb.set_label(f"log2 (clipped at ±{lim:g})", fontsize=8.5, color=MUTED_INK)
    cb.ax.tick_params(labelsize=8, colors=MUTED_INK)
    cb.outline.set_visible(False)
    fig.text(0.30, bot - 0.72 / fig.get_figheight(),
             "●  q < 0.05 within the question's family      ○  raw p < 0.05, not a hit", fontsize=8.5,
             color=TEXT_INK, va="center")
    fig.suptitle(title, x=0.01, y=1 - 0.1 / fig.get_figheight(), ha="left", fontsize=11.5, color=TEXT_INK)
    fig.text(0.01, 1 - 0.42 / fig.get_figheight(),
             "Left: model-fitted mean of each group, centred per taxon.   Right: the report's twelve comparisons.",
             ha="left", fontsize=9, color=MUTED_INK)
    save(fig, name)
    return order


# ---------------------------------------------------------------- forest
def forest(order, key, question, comps, name):
    n = len(order)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 1.5 + 0.23 * n), sharey=True)
    lo_all, hi_all = [], []
    for ax, (lab, row, st, hue) in zip(axes, comps):
        d = res[(res.row == row) & (res.stratum == st)].set_index("taxon").loc[order]
        y = np.arange(n)
        ax.axvline(0, color=RULE_INK, lw=1, zorder=0)
        for i, (_, r) in enumerate(d.iterrows()):
            hit, nom = r.q < ALPHA, r.p < ALPHA
            c = HUE[hue] if (hit or nom) else NULL_INK
            ax.plot([r.ci_lo, r.ci_hi], [i, i], color=c, lw=1.6 if (hit or nom) else 1, solid_capstyle="round")
            if hit:
                ax.plot(r.log2FC, i, "o", ms=7, color=c, mec="white", mew=1.2, zorder=3)
            elif nom:
                ax.plot(r.log2FC, i, "o", ms=6.5, mfc="white", mec=c, mew=1.6, zorder=3)
            else:
                ax.plot(r.log2FC, i, "o", ms=4, color=c, zorder=3)
        lo_all += list(d.ci_lo); hi_all += list(d.ci_hi)
        ax.set_title(lab, fontsize=10, color=HUE[hue], loc="left", fontweight="bold")
        ax.set_yticks(y, [label(t) for t in order], fontsize=8.5)
        ax.grid(axis="x", color=RULE_INK, lw=0.5, alpha=0.6)
        ax.tick_params(length=0, colors=TEXT_INK, labelsize=8.5)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(RULE_INK)
    lim = np.ceil(max(abs(np.percentile(lo_all, 1)), abs(np.percentile(hi_all, 99))))
    for ax, (lab, row, st, hue) in zip(axes, comps):
        ax.set_xlim(-lim, lim)
        dd = res[(res.row == row) & (res.stratum == st)].set_index("taxon")
        for i, t in enumerate(order):  # arrows where an interval runs off the axis
            d = dd.loc[t]
            if d.ci_lo < -lim:
                ax.annotate("", (-lim, i), (-lim + 0.6, i), arrowprops=dict(arrowstyle="->", color=NULL_INK, lw=1))
            if d.ci_hi > lim:
                ax.annotate("", (lim, i), (lim - 0.6, i), arrowprops=dict(arrowstyle="->", color=NULL_INK, lw=1))
    axes[0].invert_yaxis()
    for tl in axes[0].get_yticklabels():
        tl.set_fontstyle("italic")
    fig.supxlabel(X_LABEL[key], fontsize=9.5, color=MUTED_INK)
    fig.suptitle(question, x=0.01, ha="left", fontsize=11.5, color=TEXT_INK)
    bh = "BH within each panel" if key == "radiation" else "BH within these three"
    fig.text(0.99, 0.995, f"●  q < 0.05 ({bh})   ○  raw p < 0.05   bars: 95% interval",
             ha="right", va="top", fontsize=8.5, color=MUTED_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    save(fig, name)


# ---------------------------------------------------------------- presence
def presence_heatmap(genus_only, name):
    """Dot map: dot size = mice in the group carrying the taxon, shade = the group's mean share of reads."""
    cnt = pd.read_csv(os.path.join(TAB, f"collapsed_counts_Greengenes2_Day17.tsv"), sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(STUDY, "metadata.tsv"), sep="\t", index_col=0).loc[cnt.columns]
    depth = cnt.sum()  # collapsed_counts keeps every taxon, so these are the full depths
    cnt = cnt[(cnt > 0).sum(axis=1) >= 6]
    pres = pd.DataFrame({g: (cnt.loc[:, meta.Group == g] > 0).sum(axis=1) for g in GROUPS})
    # Mean over all four mice of the group, a mouse without the taxon counting as 0.
    rel = cnt / depth * 100
    mean = pd.DataFrame({g: rel.loc[:, meta.Group == g].mean(axis=1) for g in GROUPS})
    # One shade scale for both views, set by all 79 taxa, so a shade means the same in each.
    shown = mean.values[pres.values > 0]
    norm = LogNorm(vmin=shown.min(), vmax=shown.max())
    if genus_only:
        pres = pres[[genus_or_above(t) for t in pres.index]]
    Z = linkage(pres.values, method="average", metric="euclidean")
    order = [pres.index[i] for i in leaves_list(Z)]
    pres, mean = pres.loc[order], mean.loc[order]
    n = len(order)
    fig = plt.figure(figsize=(9.5, 2.3 + 0.2 * n))
    top = 1 - 1.3 / fig.get_figheight()
    bot = 1.0 / fig.get_figheight()
    ax_d = fig.add_axes([0.01, bot, 0.08, top - bot])
    ax = fig.add_axes([0.10, bot, 0.30, top - bot])
    dendrogram(Z, ax=ax_d, orientation="left", no_labels=True, color_threshold=0, above_threshold_color=MUTED_INK)
    ax_d.invert_yaxis()
    ax_d.axis("off")
    # Sequential blue from step 200 up, log scale: the palest dot must stay visibly blue.
    cmap = LinearSegmentedColormap.from_list("seq", ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
    area = lambda k: 150 * k / 4  # dot area proportional to mice, 4 mice fills the cell
    for i in range(n):
        ax.axhline(i, color="#eeeeeb", lw=0.5, zorder=0)
    yy, xx = np.nonzero(pres.values)
    ax.scatter(xx, yy, s=area(pres.values[yy, xx]), c=mean.values[yy, xx], cmap=cmap, norm=norm,
               edgecolors="white", linewidths=0.4, zorder=2)
    for x in (1.5, 3.5):
        ax.axvline(x, color=TEXT_INK, lw=1.2)
    ax.set_xlim(-0.5, len(GROUPS) - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(range(len(GROUPS)), GROUPS, rotation=45, ha="left", fontsize=8.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", pad=16)
    for i, g in enumerate(GROUPS):
        ax.add_patch(plt.Rectangle((i - 0.5, -1.25), 1, 0.6, color=GROUP_COLOR[g], clip_on=False))
    # No per-taxon marks: a group without the taxon shows as a missing dot.
    ax.set_yticks(range(n), [label(t) for t in order], fontsize=8.5, color=TEXT_INK)
    ax.yaxis.tick_right()
    for tl in ax.get_yticklabels():
        tl.set_fontstyle("italic")
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)

    # Shade key, ends labelled with the actual extremes so neither reads as a cap.
    h = fig.get_figheight()
    cax = fig.add_axes([0.10, bot - 0.40 / h, 0.30, 0.10 / h])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation="horizontal")
    lo, hi = norm.vmin, norm.vmax
    ticks = [t for t in (0.001, 0.01, 0.1, 1, 10)
             if lo < t < hi and abs(np.log10(t / lo)) > 0.6 and abs(np.log10(t / hi)) > 0.6]
    cb.set_ticks([lo] + ticks + [hi], labels=[f"{lo:.1g}%"] + [f"{t:g}%" for t in ticks] + [f"{hi:.2g}%"])
    cb.ax.minorticks_off()
    cb.set_label("group mean share of reads (log scale; absent mice count as 0)", fontsize=8.5, color=MUTED_INK)
    cb.ax.tick_params(labelsize=8, colors=MUTED_INK)
    cb.outline.set_visible(False)
    # Size key.
    kax = fig.add_axes([0.46, bot - 0.52 / h, 0.22, 0.30 / h])
    kax.scatter(range(4), [0] * 4, s=[area(k) for k in range(1, 5)], color=MUTED_INK, edgecolors="white", lw=0.4)
    for k in range(4):
        kax.text(k, -0.9, str(k + 1), ha="center", va="top", fontsize=8, color=MUTED_INK)
    kax.text(1.5, -2.1, "mice with the taxon, of 4 (no dot: none)", ha="center", va="top", fontsize=8.5, color=MUTED_INK)
    kax.set_xlim(-0.6, 3.6)
    kax.set_ylim(-2.5, 0.8)
    kax.axis("off")
    scope = "taxa named no further than genus" if genus_only else "taxa seen in at least 6 of 24 mice"
    fig.suptitle(f"Presence and abundance by group: the {n} {scope}", x=0.01, y=1 - 0.2 / h, ha="left",
                 fontsize=11.5, color=TEXT_INK)
    save(fig, name)


if __name__ == "__main__":
    all_taxa = list(gm.index)
    for suffix, taxa, what in (("", all_taxa, "all 46 tested taxa"),
                               ("_gplus", [t for t in all_taxa if genus_or_above(t)], "taxa named no further than genus")):
        order = abundance_heatmap(taxa, f"{STEM}{suffix}_heatmap",
                                  f"ALDEx3 abundance by group and by comparison, {len(taxa) if suffix else ''} {what}".replace(",  ", ", "))
        for key, question, comps in QUESTIONS:
            forest(order, key, question, comps, f"{STEM}{suffix}_forest_{key}")
        presence_heatmap(bool(suffix), f"presence_factorial_{TAG}{suffix}_heatmap")
