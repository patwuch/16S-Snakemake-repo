import os
# Send this job's stdout/stderr, subprocesses included, to its Snakemake log.
_log = os.open(snakemake.log[0], os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
os.dup2(_log, 1)
os.dup2(_log, 2)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from qiime2 import Artifact, Metadata
from skbio.stats.ordination import OrdinationResults

# ──────────────────────────────────────────────
# Snakemake I/O
# ──────────────────────────────────────────────
pcoa_path     = snakemake.input.pcoa
metadata_path = snakemake.input.metadata
color_palette = snakemake.params.color_palette
group_order   = snakemake.params.group_order
database      = snakemake.params.db
metric        = snakemake.params.metric
factor        = snakemake.params.group_by

output_svg      = snakemake.output.plot_svg
output_png      = snakemake.output.plot_png
output_sentinel = snakemake.output.sentinel

# ──────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────
METRIC_LABELS = {
    "jaccard":            "Jaccard",
    "braycurtis":         "Bray-Curtis",
    "weighted_unifrac":   "Weighted UniFrac",
    "unweighted_unifrac": "Unweighted UniFrac",
}
metric_label = METRIC_LABELS[metric]

# UniFrac for Silva138 is computed on a tree built by placing ASVs (via SEPP) onto
# QIIME2's SILVA 128 reference — no official SILVA 138 SEPP package exists.
# Taxonomy assignment for Silva138 is still the real SILVA 138 release; only the
# tree used for these two phylogenetic metrics is the older 128 backbone.
SILVA138_SEPP_NOTE = (
    "Note: UniFrac metrics for Silva138 use a phylogenetic tree placed via SEPP "
    "against the SILVA 128 reference (no official SILVA 138 SEPP reference "
    "exists); taxonomy assignment is still SILVA 138."
)


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    pcoa_res = Artifact.load(pcoa_path).view(OrdinationResults)

    metadata = Metadata.load(metadata_path).to_dataframe()
    for col in metadata.select_dtypes(include="object"):
        metadata[col] = metadata[col].str.strip()

    if factor not in metadata.columns:
        raise ValueError(
            f"Factor '{factor}' not found in metadata. "
            f"Available columns: {list(metadata.columns)}"
        )

    categories = group_order
    category_colors = {g: color_palette[g] for g in categories}

    coords = pcoa_res.samples.rename(columns={0: "PC1", 1: "PC2"})
    df = coords.merge(metadata, left_index=True, right_index=True)

    if df.empty:
        print(f"Warning: No overlapping indices for {metric_label}. Check your metadata IDs.")

    fig, ax = plt.subplots(figsize=(7, 6))

    sns.scatterplot(
        x="PC1",
        y="PC2",
        hue=factor,
        hue_order=categories,
        data=df,
        palette=category_colors,
        s=95,
        alpha=0.9,
        edgecolor="black",
        linewidth=0.6,
        zorder=3,
        ax=ax,
    )

    # ±1 SD ellipse along each group's own covariance axes. Visual grouping
    # only, not a confidence region.
    for cat in categories:
        cat_df = df[df[factor] == cat]
        if len(cat_df) < 3:
            continue

        vals, vecs = np.linalg.eigh(np.cov(cat_df["PC1"], cat_df["PC2"]))
        theta = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
        w, h  = 2 * np.sqrt(np.maximum(vals[::-1], 0))

        ellipse = Ellipse(
            (cat_df["PC1"].mean(), cat_df["PC2"].mean()),
            width=w,
            height=h,
            angle=theta,
            edgecolor=category_colors[cat],
            facecolor=category_colors[cat],
            alpha=0.10,
            linestyle="--",
            lw=1.0,
            zorder=1,
        )
        ax.add_patch(ellipse)

    ax.set_title(metric_label)
    ve = pcoa_res.proportion_explained.iloc[:2] * 100
    ax.set_xlabel(f"PC1 ({ve.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({ve.iloc[1]:.1f}%)")

    fig.suptitle(f"PCoA on Samples by Factor: {factor}", fontsize=14, y=0.98)
    # Small kicker line above the title, separable by cropping, so the reference
    # database used doesn't have to be inferred from the file name alone.
    fig.text(0.5, 1.05, database, ha="center", va="bottom", fontsize=9, color="gray")
    fig.tight_layout()

    if database == "Silva138" and metric in ("weighted_unifrac", "unweighted_unifrac"):
        fig.text(
            0.5, -0.02, SILVA138_SEPP_NOTE,
            ha="center", va="top", fontsize=7, style="italic", color="dimgray",
            wrap=True,
        )

    os.makedirs(os.path.dirname(output_svg), exist_ok=True)
    fig.savefig(output_svg, format="svg", dpi=300, bbox_inches="tight")
    fig.savefig(output_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(output_sentinel, "w") as fh:
        fh.write("done\n")


main()
