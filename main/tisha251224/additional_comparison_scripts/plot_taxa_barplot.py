#!/usr/bin/env python3
"""Taxa barplots (per-sample and per-group mean) for one grouping of a study's samples.

Copied from reza250909/additional_comparison_scripts/plot_taxa_barplot.py, itself from
ahmad260430/additional_comparison_scripts/plot_taxa_barplot.py (which keeps
"Other" grey instead of wrapping onto the first taxon's colour, unlike the pipeline's
plot_taxa_barplot_levels.py); --title replaces ahmad's per-arm "--arm".

Standalone counterpart to the shared pipeline's scripts/plot_taxa_barplot_levels.py.
Reuses the full study-seqs.biom table and the full per-db taxonomy TSV as-is;
restricting to the subset happens via the index intersection with the subset
metadata TSV. Must run inside the qiime2 container (needs biom + pandas).
"""
import argparse
import os
import re

import biom
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TAXA_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

CB_COLORS_20 = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
    "#AA3377", "#BBBBBB", "#000000", "#88CCEE", "#CC6677",
    "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99",
    "#999933", "#882255", "#661100", "#6699CC", "#888888",
]


OTHER_COLOR = "#E0E0E0"


def _colors(columns):
    # "Other" sits outside the palette, so with 20 named taxa it no longer wraps
    # round onto the first taxon's colour.
    return [OTHER_COLOR if c == "Other" else CB_COLORS_20[i % len(CB_COLORS_20)]
            for i, c in enumerate(columns)]


def parse_colors(pairs):
    return dict(p.split("=", 1) for p in pairs)


def _natural_sort_key(s):
    return s.map(lambda v: re.sub(r'(\d+)', lambda m: m.group().zfill(10), str(v)))


def _is_named(rank):
    """False for '', a bare prefix like 's__', or an 'Unassigned (...)' label."""
    rank = str(rank).strip()
    return not re.fullmatch(r"([a-z]__)?", rank) and not rank.lower().startswith("unassigned")


def _is_assigned(val):
    if val is None or pd.isna(val):
        return False
    return _is_named(val)


def _relabel_unassigned(row, level_index):
    parts = row["Taxon"].split(";")
    label = parts[level_index].strip() if level_index < len(parts) else ""
    if not _is_named(label):
        for higher in reversed(parts[:level_index]):
            higher = higher.strip()
            if _is_named(higher):
                return f"Unassigned ({higher})"
        return "Unassigned (All)"
    return label


def _display_label(label):
    """'f__Muribaculaceae|g__Duncaniella' -> 'Muribaculaceae | Duncaniella'.

    A species label drops its genus half when the species name repeats it.
    """
    parts = [re.sub(r"^[a-z]__", "", p.strip()) for p in str(label).split("|")]
    if len(parts) == 2 and parts[1].startswith(parts[0] + " "):
        return parts[1]
    return re.sub(r"\b[a-z]__", "", " | ".join(parts))


def _dynamic_tick_fontsize(n_samples, min_size=4, max_size=10, scale=120):
    return max(min_size, min(max_size, scale // n_samples))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-biom", required=True)
    parser.add_argument("--taxonomy-tsv", required=True)
    parser.add_argument("--metadata-tsv", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--level", required=True, choices=TAXA_LEVELS)
    parser.add_argument("--factor", required=True)
    parser.add_argument("--group-order", required=True, nargs="+")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--title", required=True, help="figure title")
    parser.add_argument("--rank-metadata-tsv", nargs="+",
                        help="choose and order the top taxa over these samples instead of this "
                             "arm's alone, so figures for different arms share taxa and colours")
    parser.add_argument("--output-samples-svg", required=True)
    parser.add_argument("--output-samples-png", required=True)
    parser.add_argument("--output-groups-svg", required=True)
    parser.add_argument("--output-groups-png", required=True)
    args = parser.parse_args()

    for out in (args.output_samples_svg, args.output_samples_png,
                args.output_groups_svg, args.output_groups_png):
        os.makedirs(os.path.dirname(out), exist_ok=True)

    level = args.level
    factor = args.factor
    group_order = args.group_order
    db = args.db
    top_n = args.top_n

    table = biom.load_table(args.table_biom)
    df = table.to_dataframe(dense=True).T
    df.index = df.index.astype(str).str.strip()

    taxonomy = pd.read_csv(args.taxonomy_tsv, sep="\t", index_col=0)
    df = df[df.columns.intersection(taxonomy.index)]
    taxonomy = taxonomy.loc[df.columns]

    metadata = pd.read_csv(args.metadata_tsv, sep="\t", index_col=0)
    metadata.index = metadata.index.astype(str).str.strip()

    rank_samples = metadata.index
    if args.rank_metadata_tsv:
        rank_samples = pd.Index([])
        for path in args.rank_metadata_tsv:
            rank_samples = rank_samples.union(pd.read_csv(path, sep="\t", index_col=0).index.astype(str).str.strip())
    df_rank = df.loc[df.index.intersection(rank_samples)]

    common_samples = df.index.intersection(metadata.index)
    df = df.loc[common_samples]
    metadata = metadata.loc[common_samples]

    tax_split = taxonomy["Taxon"].str.split(";", expand=True)
    for i, lvl in enumerate(TAXA_LEVELS):
        taxonomy[lvl] = tax_split[i].str.strip() if i in tax_split.columns else ""

    level_index = TAXA_LEVELS.index(level)
    if level == "Kingdom":
        taxonomy[level] = taxonomy["Taxon"].str.split(";").str[level_index].str.strip().fillna("Unassigned (All)")
    else:
        taxonomy[level] = taxonomy.apply(_relabel_unassigned, axis=1, level_index=level_index)

    if level == "Genus":
        taxonomy[level] = taxonomy.apply(
            lambda r: f"{r['Family']}|{r[level]}"
            if _is_assigned(r.get("Family")) and _is_assigned(r.get(level)) else r[level],
            axis=1,
        )
    elif level == "Species":
        taxonomy[level] = taxonomy.apply(
            lambda r: f"{r['Genus']}|{r[level]}"
            if _is_assigned(r.get("Genus")) and _is_assigned(r.get(level)) else r[level],
            axis=1,
        )

    metadata["_group_rank"] = metadata[factor].map({g: i for i, g in enumerate(group_order)})
    df = df.loc[metadata.sort_values(by=["_group_rank"], kind="stable", key=_natural_sort_key).index]
    metadata = metadata.loc[df.index]

    df_tax = df.T.groupby(taxonomy[level]).sum().T
    df_tax_norm = df_tax.div(df_tax.sum(axis=1), axis=0)

    df_rank_tax = df_rank.T.groupby(taxonomy[level]).sum().T
    df_rank_norm = df_rank_tax.div(df_rank_tax.sum(axis=1), axis=0)
    top_taxa = df_rank_norm.mean().sort_values(ascending=False).head(top_n).index
    top_taxa = top_taxa.intersection(df_tax_norm.columns, sort=False)
    df_plot = df_tax_norm[top_taxa].copy()
    df_plot["Other"] = df_tax_norm.drop(columns=top_taxa).sum(axis=1)

    n_samples = len(df_plot)
    tick_fontsize = _dynamic_tick_fontsize(n_samples)
    colors = _colors(df_plot.columns)

    fig, ax = plt.subplots(figsize=(12, 6))
    bottom = np.zeros(n_samples)
    for i, col in enumerate(df_plot.columns):
        ax.bar(df_plot.index, df_plot[col], bottom=bottom, color=colors[i], label=col, width=0.8, linewidth=0)
        bottom += df_plot[col].values

    group_ordered = metadata.loc[df_plot.index, factor]
    y_max = float(np.nan_to_num(bottom.max(), nan=1.0, posinf=1.0))
    ax.set_ylim(0, y_max * 1.08)

    for group in group_order:
        idxs = np.where(group_ordered == group)[0]
        if len(idxs) == 0:
            continue
        if idxs[0] > 0:
            ax.axvline(idxs[0] - 0.5, color="gray", linestyle="--", linewidth=0.7, alpha=0.5)
        ax.text((idxs[0] + idxs[-1]) / 2, y_max * 1.01, str(group), ha="center", va="bottom",
                fontsize=10, fontweight="bold", clip_on=False)

    ax.set_ylabel("Relative abundance")
    ax.set_xlabel("Samples")
    ax.set_title(args.title)
    ax.annotate(f"{db} · {level}", xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 18),
                textcoords="offset points", ha="center", va="bottom", fontsize=9, color="gray")

    plt.xticks(ticks=np.arange(n_samples), labels=df_plot.index, rotation=90, fontsize=tick_fontsize)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], [_display_label(l) for l in labels[::-1]], bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8, frameon=False)

    plt.tight_layout()
    plt.savefig(args.output_samples_svg, format="svg", bbox_inches="tight")
    plt.savefig(args.output_samples_png, format="png", dpi=300, bbox_inches="tight")
    plt.close()

    df_plot["_group"] = metadata.loc[df_plot.index, factor]
    df_grouped = df_plot.groupby("_group")[list(top_taxa) + ["Other"]].mean().loc[group_order]

    n_groups = len(df_grouped)
    colors = _colors(df_grouped.columns)

    fig, ax = plt.subplots(figsize=(max(6, n_groups * 1.5), 6))
    bottom = np.zeros(n_groups)
    for i, col in enumerate(df_grouped.columns):
        ax.bar(df_grouped.index, df_grouped[col], bottom=bottom, color=colors[i], label=col, width=0.6, linewidth=0)
        bottom += df_grouped[col].values

    y_max = float(np.nan_to_num(bottom.max(), nan=1.0, posinf=1.0))
    ax.set_ylim(0, y_max * 1.05)
    ax.set_ylabel("Mean relative abundance")
    ax.set_xlabel(factor)
    ax.set_title(args.title)
    ax.annotate(f"{db} · {level}", xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 18),
                textcoords="offset points", ha="center", va="bottom", fontsize=9, color="gray")

    plt.xticks(ticks=np.arange(n_groups), labels=df_grouped.index, rotation=45, ha="right", fontsize=10)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], [_display_label(l) for l in labels[::-1]], bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8, frameon=False)

    plt.tight_layout()
    plt.savefig(args.output_groups_svg, format="svg", bbox_inches="tight")
    plt.savefig(args.output_groups_png, format="png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {args.output_samples_png}, {args.output_groups_png}")


if __name__ == "__main__":
    main()
