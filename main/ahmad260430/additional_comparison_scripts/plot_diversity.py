#!/usr/bin/env python3
"""Alpha/beta diversity box- and PCoA-plots for a Day17-vs-Day25 subset.

Standalone counterpart to the shared pipeline's scripts/plot_alpha_diversity_metric.py
and scripts/plot_beta_diversity_metric.py. Reuses the already-computed
core-metrics vectors/PCoAs for the *full* sample set — no recomputation needed,
since restricting to the subset happens by joining against the subset
metadata (alpha: left-join, unmatched samples get NaN and are dropped by the
categorical plot; beta: inner-merge). Must run inside the qiime2 container.
"""
import argparse
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from qiime2 import Artifact, Metadata
from skbio.stats.ordination import OrdinationResults

ALPHA_METRIC_LABELS = {
    "shannon": "Shannon", "evenness": "Pielou's Evenness",
    "simpson": "Simpson", "chao1": "Chao1", "faith_pd": "Faith's PD",
}
BETA_METRIC_LABELS = {
    "jaccard": "Jaccard", "braycurtis": "Bray-Curtis",
    "weighted_unifrac": "Weighted UniFrac", "unweighted_unifrac": "Unweighted UniFrac",
}


def parse_colors(pairs):
    return dict(p.split("=", 1) for p in pairs)


def mouse_id(sid):
    """E1_1.00 and 1.00colon are the same mouse at two timepoints."""
    m = re.match(r"^E\d+_(\d+)\.(\d+)$", sid) or re.match(r"^(\d+)\.(\d+)colon$", sid)
    return f"{m.group(1)}_{m.group(2)}" if m else None


def plot_alpha(args):
    metric_label = ALPHA_METRIC_LABELS[args.metric]
    vector_df = Artifact.load(args.vector).view(pd.Series).rename(metric_label).to_frame()
    metadata = Metadata.load(args.metadata).to_dataframe()
    merged = vector_df.join(metadata)

    colors = parse_colors(args.color)
    group_order = args.group_order
    palette = [colors[g] for g in group_order]

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(max(4, len(group_order) * 1.2), 5))
    sns.boxplot(data=merged, x=args.factor, y=metric_label, order=group_order, palette=palette, ax=ax)

    # One line per mouse across the timepoints, so the pairing is visible at n = 4.
    # Sample ids that do not carry a mouse fall back to the unpaired scatter.
    points = merged.dropna(subset=[metric_label, args.factor]).copy()
    points["MouseID"] = [mouse_id(s) for s in points.index]
    xpos = {g: i for i, g in enumerate(group_order)}
    if points["MouseID"].notna().any():
        for _, rows in points.dropna(subset=["MouseID"]).groupby("MouseID"):
            pts = sorted((xpos[g], v) for g, v in zip(rows[args.factor], rows[metric_label]) if g in xpos)
            if not pts:
                continue
            xs, ys = zip(*pts)
            if len(pts) > 1:
                ax.plot(xs, ys, color="0.25", lw=1.1, alpha=0.75, zorder=3)
            ax.plot(xs, ys, ls="none", marker="o", ms=4.5, zorder=4,
                    mfc="white", mec="0.25", mew=1.1)
    else:
        sns.stripplot(data=merged, x=args.factor, y=metric_label, order=group_order,
                      color="black", alpha=0.5, ax=ax)

    tick_labels = ax.get_xticklabels()
    if tick_labels:
        ax.set_xticklabels([t.get_text() for t in tick_labels], rotation=25, ha="right", rotation_mode="anchor")
    ax.set_xlabel(args.factor)
    ax.set_ylabel(metric_label)
    ax.set_title(f"{metric_label} by {args.factor}")
    if args.db:
        ax.annotate(args.db, xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 18),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=9, color="gray")

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output_svg), exist_ok=True)
    fig.savefig(args.output_svg, format="svg", bbox_inches="tight")
    fig.savefig(args.output_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {args.output_svg}")


def plot_beta(args):
    metric_label = BETA_METRIC_LABELS[args.metric]
    pcoa_res = Artifact.load(args.pcoa).view(OrdinationResults)
    metadata = Metadata.load(args.metadata).to_dataframe()
    for col in metadata.select_dtypes(include="object"):
        metadata[col] = metadata[col].str.strip()

    if args.factor not in metadata.columns:
        raise ValueError(f"Factor '{args.factor}' not found in metadata. Available: {list(metadata.columns)}")

    colors = parse_colors(args.color)
    categories = args.group_order
    category_colors = {g: colors[g] for g in categories}

    coords = pcoa_res.samples.rename(columns={0: "PC1", 1: "PC2"})
    df = coords.merge(metadata, left_index=True, right_index=True)
    if df.empty:
        print(f"Warning: No overlapping indices for {metric_label}. Check your metadata IDs.")

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(x="PC1", y="PC2", hue=args.factor, hue_order=categories, data=df,
                     palette=category_colors, s=100, alpha=0.8, ax=ax)

    for cat in categories:
        cat_df = df[df[args.factor] == cat]
        if len(cat_df) < 2:
            continue
        if cat_df["PC1"].std() == 0 or cat_df["PC2"].std() == 0:
            continue
        cx, cy = cat_df["PC1"].mean(), cat_df["PC2"].mean()
        w, h = cat_df["PC1"].std() * 2, cat_df["PC2"].std() * 2
        ellipse = Ellipse((cx, cy), width=w, height=h, edgecolor=category_colors[cat],
                           facecolor=category_colors[cat], alpha=0.15, lw=2)
        ax.add_patch(ellipse)
        ax.scatter(cx, cy, marker="D", s=150, color=category_colors[cat])

    ax.set_title(metric_label)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    fig.suptitle(f"PCoA on Samples by Factor: {args.factor}", fontsize=14, y=0.98)
    if args.db:
        fig.text(0.5, 1.05, args.db, ha="center", va="bottom", fontsize=9, color="gray")
    fig.tight_layout()

    os.makedirs(os.path.dirname(args.output_svg), exist_ok=True)
    fig.savefig(args.output_svg, format="svg", dpi=300, bbox_inches="tight")
    fig.savefig(args.output_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {args.output_svg}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)

    a = sub.add_parser("alpha")
    a.add_argument("--vector", required=True)
    a.add_argument("--metadata", required=True)
    a.add_argument("--metric", required=True, choices=sorted(ALPHA_METRIC_LABELS))
    a.add_argument("--factor", required=True)
    a.add_argument("--db", default="")
    a.add_argument("--group-order", required=True, nargs="+")
    a.add_argument("--color", action="append", required=True, metavar="GROUP=HEX")
    a.add_argument("--output-svg", required=True)
    a.add_argument("--output-png", required=True)
    a.set_defaults(func=plot_alpha)

    b = sub.add_parser("beta")
    b.add_argument("--pcoa", required=True)
    b.add_argument("--metadata", required=True)
    b.add_argument("--metric", required=True, choices=sorted(BETA_METRIC_LABELS))
    b.add_argument("--factor", required=True)
    b.add_argument("--db", default="")
    b.add_argument("--group-order", required=True, nargs="+")
    b.add_argument("--color", action="append", required=True, metavar="GROUP=HEX")
    b.add_argument("--output-svg", required=True)
    b.add_argument("--output-png", required=True)
    b.set_defaults(func=plot_beta)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
