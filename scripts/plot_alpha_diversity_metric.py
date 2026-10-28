import os
# Send this job's stdout/stderr, subprocesses included, to its Snakemake log.
_log = os.open(snakemake.log[0], os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
os.dup2(_log, 1)
os.dup2(_log, 2)

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from qiime2 import Artifact, Metadata

# ──────────────────────────────────────────────
# Snakemake I/O
# ──────────────────────────────────────────────
vector_path   = snakemake.input.vector
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
    "shannon":  "Shannon",
    "evenness": "Pielou's Evenness",
    "simpson":  "Simpson",
    "chao1":    "Chao1",
    "faith_pd": "Faith's PD",
}
metric_label = METRIC_LABELS[metric]

# Faith PD for Silva138 is computed on a tree built by placing ASVs (via SEPP)
# onto QIIME2's SILVA 128 reference — no official SILVA 138 SEPP package exists.
# Taxonomy assignment for Silva138 is still the real SILVA 138 release; only the
# tree used for this one phylogenetic metric is the older 128 backbone.
SILVA138_SEPP_NOTE = (
    "Note: Faith PD for Silva138 uses a phylogenetic tree placed via SEPP against "
    "the SILVA 128 reference (no official SILVA 138 SEPP reference exists); "
    "taxonomy assignment is still SILVA 138."
)


# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────
def load_vector(path, label):
    return Artifact.load(path).view(pd.Series).rename(label).to_frame()


def load_metadata(path):
    return Metadata.load(path).to_dataframe()


def plot_metric(df, factor_col, metric_label, outfile_svg, outfile_png, color_palette, group_order):
    """Box + strip plot for a single alpha-diversity metric, saved as a standalone figure."""
    palette = [color_palette[g] for g in group_order]

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(max(4, len(group_order) * 1.2), 5))

    sns.boxplot(
        data=df,
        x=factor_col,
        y=metric_label,
        order=group_order,
        palette=palette,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x=factor_col,
        y=metric_label,
        order=group_order,
        color="black",
        alpha=0.5,
        ax=ax,
    )

    tick_labels = ax.get_xticklabels()
    if tick_labels:
        ax.set_xticklabels(
            [t.get_text() for t in tick_labels],
            rotation=25,
            ha="right",
            rotation_mode="anchor",
        )
    ax.set_xlabel(factor_col)
    ax.set_ylabel(metric_label)
    ax.set_title(f"{metric_label} by {factor_col}")
    # Small kicker line above the title, separable by cropping, so the reference
    # database used doesn't have to be inferred from the file name alone.
    ax.annotate(
        database,
        xy=(0.5, 1.0), xycoords="axes fraction",
        xytext=(0, 18), textcoords="offset points",
        ha="center", va="bottom", fontsize=9, color="gray",
    )

    if database == "Silva138" and metric == "faith_pd":
        fig.text(
            0.5, -0.02, SILVA138_SEPP_NOTE,
            ha="center", va="top", fontsize=7, style="italic", color="dimgray",
            wrap=True,
        )

    fig.tight_layout()
    fig.savefig(outfile_svg, format="svg", bbox_inches="tight")
    fig.savefig(outfile_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    vector_df = load_vector(vector_path, metric_label)
    metadata  = load_metadata(metadata_path)
    merged    = vector_df.join(metadata)

    os.makedirs(os.path.dirname(output_svg), exist_ok=True)

    plot_metric(merged, factor, metric_label, output_svg, output_png, color_palette, group_order)
    print(f"Saved: {output_svg}")

    with open(output_sentinel, "w") as fh:
        fh.write("done\n")


main()
