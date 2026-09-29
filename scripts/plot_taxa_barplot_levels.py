import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import biom

# ================================
# INPUTS FROM SNAKEMAKE
# ================================
table_biom       = snakemake.input.table_biom
taxonomy_tsv     = snakemake.input.taxonomy_tsv
metadata_tsv     = snakemake.params.metadata_tsv
db               = snakemake.params.database
top_n            = snakemake.params.top_n_taxa_shown_on_barplot
level            = snakemake.params.taxa_level
factor           = snakemake.params.group_by
group_order      = snakemake.params.group_order
output_path_samples     = snakemake.output.plot_samples
output_path_groups      = snakemake.output.plot_groups
output_path_samples_png = snakemake.output.plot_samples_png
output_path_groups_png  = snakemake.output.plot_groups_png
dropped_samples  = set(map(str, snakemake.params.get("dropped_sampleid", [])))

# ================================
# CONSTANTS
# ================================
TAXA_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

CB_COLORS_20 = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
    "#AA3377", "#BBBBBB", "#000000", "#88CCEE", "#CC6677",
    "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99",
    "#999933", "#882255", "#661100", "#6699CC", "#888888",
]
OTHER_COLOR = "#E0E0E0"

# ================================
# HELPER FUNCTIONS
# ================================

def _natural_sort_key(s):
    """Zero-pad embedded integers so '10' sorts after '5'."""
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
    """Return a human-readable label, falling back to higher-rank info when unassigned."""
    parts = row["Taxon"].split(";")
    label = parts[level_index].strip() if level_index < len(parts) else ""
    if not _is_named(label):
        for higher in reversed(parts[:level_index]):
            higher = higher.strip()
            if _is_named(higher):
                return f"Unassigned ({higher})"
        return "Unassigned (All)"
    return label


def _colors(columns):
    # "Other" sits outside the palette, so with 20 named taxa it no longer wraps
    # round onto the first taxon's colour.
    return [OTHER_COLOR if c == "Other" else CB_COLORS_20[i % len(CB_COLORS_20)]
            for i, c in enumerate(columns)]


def _dynamic_tick_fontsize(n_samples, min_size=4, max_size=10, scale=120):
    """Return a font size that shrinks as sample count grows."""
    return max(min_size, min(max_size, scale // n_samples))


# ================================
# LOAD & ALIGN DATA
# ================================
table = biom.load_table(table_biom)
df = table.to_dataframe(dense=True).T          # samples × features
df.index = df.index.astype(str).str.strip()

taxonomy = pd.read_csv(taxonomy_tsv, sep="\t", index_col=0)
df       = df[df.columns.intersection(taxonomy.index)]
taxonomy = taxonomy.loc[df.columns]

metadata = pd.read_csv(metadata_tsv, sep="\t", index_col=0)
metadata.index = metadata.index.astype(str).str.strip()

common_samples = df.index.intersection(metadata.index)
df       = df.loc[common_samples]
metadata = metadata.loc[common_samples]

# ================================
# DROP UNWANTED SAMPLES
# ================================
if dropped_samples:
    keep = ~df.index.isin(dropped_samples)
    df       = df.loc[keep]
    metadata = metadata.loc[df.index]

# ================================
# PARSE TAXONOMY
# ================================
tax_split = taxonomy["Taxon"].str.split(";", expand=True)
for i, lvl in enumerate(TAXA_LEVELS):
    taxonomy[lvl] = tax_split[i].str.strip() if i in tax_split.columns else ""

level_index = TAXA_LEVELS.index(level)

# Assign readable labels at the requested level
if level == "Kingdom":
    taxonomy[level] = (
        taxonomy["Taxon"].str.split(";").str[level_index].str.strip().fillna("Unassigned (All)")
    )
else:
    taxonomy[level] = taxonomy.apply(_relabel_unassigned, axis=1, level_index=level_index)


# # ============== Bacillota collapse (Optional, comment out if undesired)==============
# def _merge_bacillota(label):
#     """Collapse all GTDB Bacillota sub-phyla (Bacillota, Bacillota_A, Bacillota_B, …)
#     into a single label, preserving any rank prefix like 'p__'."""
#     m = re.match(r'^([a-zA-Z]__)?Bacillota', str(label).strip())
#     if m:
#         return f"{m.group(1) or ''}Bacillota"
#     return label


# if level == "Phylum":
#     taxonomy[level] = taxonomy[level].map(_merge_bacillota)
# # =================================================================================


# Prefix Genus / Species labels with their parent rank for disambiguation
if level == "Genus":
    taxonomy[level] = taxonomy.apply(
        lambda r: f"{r['Family']}|{r[level]}"
        if _is_assigned(r.get("Family")) and _is_assigned(r.get(level))
        else r[level],
        axis=1,
    )
elif level == "Species":
    taxonomy[level] = taxonomy.apply(
        lambda r: f"{r['Genus']}|{r[level]}"
        if _is_assigned(r.get("Genus")) and _is_assigned(r.get(level))
        else r[level],
        axis=1,
    )

# ================================
# ORDER SAMPLES
# ================================
factor_list   = [f.strip() for f in factor.split(",")]
primary_factor = factor_list[0]

# Build a sort key: first by group position in group_order, then by remaining factors
metadata["_group_rank"] = metadata[primary_factor].map(
    {g: i for i, g in enumerate(group_order)}
)
sort_cols = ["_group_rank"] + factor_list[1:]
df = df.loc[
    metadata.sort_values(by=sort_cols, kind="stable", key=_natural_sort_key).index
]
metadata = metadata.loc[df.index]

# ================================
# AGGREGATE & NORMALISE
# ================================
df_tax      = df.T.groupby(taxonomy[level]).sum().T
df_tax_norm = df_tax.div(df_tax.sum(axis=1), axis=0)

# ================================
# SELECT TOP N + "OTHER"
# ================================
top_taxa = df_tax_norm.mean().sort_values(ascending=False).head(top_n).index
df_plot  = df_tax_norm[top_taxa].copy()
df_plot["Other"] = df_tax_norm.drop(columns=top_taxa).sum(axis=1)

# ================================
# PLOT
# ================================
n_samples     = len(df_plot)
tick_fontsize = _dynamic_tick_fontsize(n_samples)
colors        = _colors(df_plot.columns)

fig, ax = plt.subplots(figsize=(12, 6))

bottom = np.zeros(n_samples)
for i, col in enumerate(df_plot.columns):
    ax.bar(
        df_plot.index,
        df_plot[col],
        bottom=bottom,
        color=colors[i],
        label=col,
        width=0.8,
        linewidth=0,
    )
    bottom += df_plot[col].values

# Group annotations — iterate in Order-derived sequence
group_ordered = metadata.loc[df_plot.index, primary_factor]
y_max = float(np.nan_to_num(bottom.max(), nan=1.0, posinf=1.0))
ax.set_ylim(0, y_max * 1.08)

for group in group_order:
    idxs = np.where(group_ordered == group)[0]
    if len(idxs) == 0:
        continue
    if idxs[0] > 0:
        ax.axvline(idxs[0] - 0.5, color="gray", linestyle="--", linewidth=0.7, alpha=0.5)
    ax.text(
        (idxs[0] + idxs[-1]) / 2,
        y_max * 1.01,
        str(group),
        ha="center", va="bottom",
        fontsize=10, fontweight="bold",
        clip_on=False,
    )

# Axes labels & title
ax.set_ylabel("Relative abundance")
ax.set_xlabel("Samples")
ax.set_title(f"Relative abundance at {level} level")
# Small kicker line above the title, separable by cropping, so the reference
# database used doesn't have to be inferred from the file name alone.
ax.annotate(
    db,
    xy=(0.5, 1.0), xycoords="axes fraction",
    xytext=(0, 18), textcoords="offset points",
    ha="center", va="bottom", fontsize=9, color="gray",
)

plt.xticks(
    ticks=np.arange(n_samples),
    labels=df_plot.index,
    rotation=90,
    fontsize=tick_fontsize,
)

handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[::-1], labels[::-1],
    bbox_to_anchor=(1.05, 1), loc="upper left",
    fontsize=8, frameon=False,
)

plt.tight_layout()
plt.savefig(output_path_samples, format="svg", bbox_inches="tight")
plt.savefig(output_path_samples_png, format="png", dpi=300, bbox_inches="tight")
plt.close()

# ================================
# AGGREGATE BY GROUP
# ================================
df_plot["_group"] = metadata.loc[df_plot.index, primary_factor]

# Mean relative abundance per group, preserving group_order
df_grouped = (
    df_plot
    .groupby("_group")[list(top_taxa) + ["Other"]]
    .mean()
    .loc[group_order]          # enforce the Order-derived sequence
)

# ================================
# PLOT GROUP MEANS
# ================================
n_groups      = len(df_grouped)
colors        = _colors(df_grouped.columns)

fig, ax = plt.subplots(figsize=(max(6, n_groups * 1.5), 6))

bottom = np.zeros(n_groups)
for i, col in enumerate(df_grouped.columns):
    ax.bar(
        df_grouped.index,
        df_grouped[col],
        bottom=bottom,
        color=colors[i],
        label=col,
        width=0.6,
        linewidth=0,
    )
    bottom += df_grouped[col].values

y_max = float(np.nan_to_num(bottom.max(), nan=1.0, posinf=1.0))
ax.set_ylim(0, y_max * 1.05)

ax.set_ylabel("Mean relative abundance")
ax.set_xlabel(primary_factor)
ax.set_title(f"Relative abundance at {level} level")
# Small kicker line above the title, separable by cropping, so the reference
# database used doesn't have to be inferred from the file name alone.
ax.annotate(
    db,
    xy=(0.5, 1.0), xycoords="axes fraction",
    xytext=(0, 18), textcoords="offset points",
    ha="center", va="bottom", fontsize=9, color="gray",
)

plt.xticks(ticks=np.arange(n_groups), labels=df_grouped.index, rotation=45, ha="right", fontsize=10)

handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[::-1], labels[::-1],
    bbox_to_anchor=(1.05, 1), loc="upper left",
    fontsize=8, frameon=False,
)

plt.tight_layout()
plt.savefig(output_path_groups, format="svg", bbox_inches="tight")
plt.savefig(output_path_groups_png, format="png", dpi=300, bbox_inches="tight")
plt.close()