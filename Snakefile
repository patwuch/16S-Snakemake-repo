import os
import re
import sys
import glob
import warnings
from pathlib import Path
import pandas as pd
from snakemake.io import directory

def log(msg):
    print(msg, file=sys.stderr, flush=True)

def reject_whitespace(value, what):
    """Reject whitespace in a value that becomes a shell-visible path segment.

    This applies to STUDY_NAME and nothing else. STUDY_DIR sits under it, and
    several rules interpolate STUDY_DIR-relative paths into unquoted shell --
    `mkdir -p {params.outdir}` in rules/NGS.smk and rules/TGS.smk, among others
    -- where a space silently splits one path into two arguments.

    Grouping COLUMN names are deliberately not restricted: snakemake
    round-trips such a name through wildcards, and no shell rule ever carries a
    grouping column in its paths.

    Group VALUES are not restricted either -- "Low tumor volume" survives end to
    end by design.
    """
    if any(c.isspace() for c in str(value)):
        raise ValueError(
            f"{what} must not contain whitespace: {value!r}. It becomes a "
            "directory and file-name segment that reaches unquoted shell "
            "commands. Use underscores instead. Note this restricts the study "
            "name only -- metadata column names and group values may contain "
            "spaces."
        )
    return value

# Suppress gRPC C-level verbosity (controls gRPC's own log system)
os.environ["GRPC_VERBOSITY"] = "NONE"
# Suppress abseil (absl) logs — this covers the E0000 instrument.cc duplicate-metric
# registration messages that GRPC_VERBOSITY alone does not silence
os.environ["GLOG_minloglevel"] = "3"
# Suppress emperor's pkg_resources deprecation warning (emperor 1.0.4 uses the
# deprecated pkg_resources API; setuptools 67.5+ warns regardless of version pin)
warnings.filterwarnings("ignore", "pkg_resources is deprecated", UserWarning)
shell.prefix('export PATH="$CONDA_PREFIX/bin:$PATH"; ')


###############################
# ==============================
# Load config parameters
# ==============================

# --- Sequence type and region validation ---
sequence_type = config.get('SEQUENCE_TYPE')
if sequence_type not in ["NGS", "TGS"]:
    raise ValueError("SEQUENCE_TYPE must be 'NGS' or 'TGS'.")

region = config.get("REGION")
if region not in ['region_V3V4', 'full_length']:
    raise ValueError(f"{region} not supported. Please choose a valid region.")

study_name = reject_whitespace(config["STUDY_NAME"], "STUDY_NAME")

# --- Analysis mode ---
ANALYSIS_MODE = config.get("ANALYSIS_MODE", "standard")
if ANALYSIS_MODE not in ["standard"]:
    raise ValueError(f"Unknown ANALYSIS_MODE: {ANALYSIS_MODE}")

# Rules cheap enough that submitting them as their own SLURM job (via
# --workflow-profile profiles/slurm) would cost more in sbatch latency than
# the rule itself takes to run.
localrules: generate_manifest

# ==============================
# Experimental design validation
# ==============================
def validate_standard_design(config):
    log("[VALIDATION] Checking experimental design (standard mode)...")
    design = config.get("DESIGN", {})

    # This pipeline fits the flattest model there is: one factor column per run,
    # as a flat categorical. Terms that adjust or condition that fit are rejected
    # here rather than quietly dropped, because a key that silently does nothing
    # is how DESIGN.INTERACTIONS sat broken and unnoticed for as long as it did.
    # Whether a column is a coefficient or a nuisance confounder is a property of
    # the study, and it changes faster than a shared schema can track.
    unsupported = sorted(k for k in ("COVARIATES", "INTERACTIONS") if design.get(k))
    if unsupported:
        raise ValueError(
            f"DESIGN.{', DESIGN.'.join(unsupported)} is set, and this pipeline does "
            "not model it. Every analysis treats one DESIGN.FACTORS column at a "
            "time as a flat categorical. Formula-based modelling (interactions, covariate adjustment, "
            "random effects) belongs in a per-study script under "
            "main/<study>/additional_comparison_scripts/, where the modelling "
            "choice sits next to the reasoning for it."
        )

    factors_raw = design.get("FACTORS") or []
    factors = [f["column"] for f in factors_raw]

    composites_raw = design.get("COMPOSITE_LABELS") or []
    composites = [c["column"] for c in composites_raw]

    # Basic sanity checks
    if not factors:
        raise ValueError("Standard mode requires at least one DESIGN.FACTOR.")

    # Check uniqueness across all columns
    all_columns = factors + composites
    if len(all_columns) != len(set(all_columns)):
        raise ValueError("Duplicate column names detected across FACTORS or COMPOSITE_LABELS.")

    # Ensure composites are not used as model factors
    for comp in composites:
        if comp in factors:
            raise ValueError(f"Composite label '{comp}' must not be used as a FACTOR.")

    log(f"✓ Factors: {factors}")
    log(f"✓ Composite labels (visualization only): {composites}")

    return {
        "factors": factors,
        "composite_labels": composites,
    }
# --- Determine grouping axes based on analysis mode ---
DESIGN_INFO = config.get("DESIGN", {})

if ANALYSIS_MODE == "standard":
    # Validate full design
    validated_info = validate_standard_design(config)
    GROUPING_AXES = validated_info.get("factors", [])

    if not GROUPING_AXES:
        raise ValueError("Standard mode requires at least one factor.")

    # Include grouping axes in DESIGN_INFO for downstream use
    validated_info["mode"] = "standard"
    validated_info["grouping_axes"] = GROUPING_AXES
    DESIGN_INFO = validated_info

else:
    raise ValueError(f"Unknown ANALYSIS_MODE: {ANALYSIS_MODE}")

log(f"[MAIN] Grouping axes for outputs: {GROUPING_AXES}")

# ==============================
# Reference DB validation
# ==============================
allowed_reference_dbs = ["Greengenes2", "Silva138"]
reference_db = config.get("REFERENCE_DB", [])
for db in reference_db:
    if db not in allowed_reference_dbs:
        raise ValueError(f"Invalid REFERENCE_DB '{db}'. Allowed: {allowed_reference_dbs}")

# Silva138's tree is built per-study via SEPP fragment-insertion (memory-heavy,
# has previously OOM-killed this host — see rules/Silva138.smk). When skipped,
# targets that depend on a per-db tree (Faith PD, UniFrac, their PCoA, their
# rarefaction depth) are built over this narrower list instead of reference_db;
# everything else (taxonomy, taxa barplots, non-phylogenetic diversity, etc.)
# keeps using reference_db unchanged.
skip_silva138_phylogeny = config.get("SILVA138_SKIP_PHYLOGENY", False)
phylogenetic_reference_db = [
    db for db in reference_db
    if not (db == "Silva138" and skip_silva138_phylogeny)
]

# ==============================
# Infer differential abundance method
# ==============================
DA_METHODS = []
if ANALYSIS_MODE == "standard":
    GROUPING_AXES = validated_info.get("factors", [])
    validated_info["mode"] = "standard"
    validated_info["grouping_axes"] = GROUPING_AXES
    DESIGN_INFO = validated_info
    # --- UPDATED DA LOGIC ---
    if DESIGN_INFO.get("factors"):
        # Automatically add LEfSe for each factor if desired
        DA_METHODS.append("LEfSe_per_factor")


# ==============================
# Per-implementation toggles (LEfSe, PICRUSt2)
# ==============================
# Each analysis is switched independently; see pipeline_overview.md Stage 5 for
# what each one produces. Both default to True (every method always ran
# unconditionally before these toggles existed) so configs that predate this
# option keep running every method unless they opt out.
#
# Two DA methods are gone, and a leftover toggle for either in an old config is
# simply ignored:
#   RUN_LEFSE_MICROBIOMEMARKER - the microbiomeMarker LEfSe reimplementation added
#     little over the Huttenhower original it reimplements.
#   RUN_ANCOMBC2 (and the ANCOMBC2_* parameters) - ANCOM-BC2's strength is
#     modelling a study's design, which a project-agnostic flat fit cannot do;
#     its calls here were easy to misread. It runs in per-study scripts instead.
run_lefse_huttenhower      = config.get("RUN_LEFSE_HUTTENHOWER", True)
run_picrust2                = config.get("RUN_PICRUST2", True)

_enabled_da_methods = [name for name, on in [
    ("lefse_Huttenhower", run_lefse_huttenhower),
] if on]
log(f"[MAIN] Differential abundance methods to be run: {', '.join(_enabled_da_methods) if _enabled_da_methods else 'none'}")
log(f"[MAIN] PICRUSt2: {'enabled' if run_picrust2 else 'disabled'}")

# ==============================
# Other parameters
# ==============================
n_threads = config.get("N_THREADS", 12)

# ==============================
# Advanced parameters 
# ==============================
cutadapt_error_rate          = config.get("CUTADAPT_ERROR_RATE", 0.1)
dada2_trim_left_f            = config.get("DADA2_TRIM_LEFT_F", 0)
dada2_trim_left_r            = config.get("DADA2_TRIM_LEFT_R", 0)
trunc_len_q_threshold        = config.get("TRUNC_LEN_Q_THRESHOLD", 20)
tgs_dada2_trunc_len          = config.get("TGS_DADA2_TRUNC_LEN", 0)
rarefy_depth_percentile      = config.get("RAREFY_DEPTH_PERCENTILE", 0)
top_n_taxa                   = config.get("TOP_N_TAXA", 20)
top_n_picrust                = config.get("TOP_N_PICRUST", 30)
lefse_lda_cutoff             = config.get("LEFSE_LDA_CUTOFF", 3.0)
lefse_kw_cutoff              = config.get("LEFSE_KW_CUTOFF", 0.05)
lefse_wilcoxon_cutoff        = config.get("LEFSE_WILCOXON_CUTOFF", 0.05)

# ==============================
# Summary
# ==============================
log("="*50)
log("[MAIN] Config loaded successfully")
log(f"[MAIN] sequence_type = '{sequence_type}'")
log(f"[MAIN] reference_db = {reference_db}")
log(f"[MAIN] study_name = '{study_name}'")
log(f"[MAIN] Analysis mode = '{ANALYSIS_MODE}'")
log(f"[MAIN] Using grouping: {DESIGN_INFO.get('grouping_column', DESIGN_INFO.get('factors', []))}")
log(f"[MAIN] About to include: rules/{sequence_type}.smk")
log("="*50)

# Global directories
SNAKEFILE_DIR = Path.cwd()  # current working directory
REPO_ROOT = SNAKEFILE_DIR
WORKFLOW_DIR = SNAKEFILE_DIR
MAIN_DIR = REPO_ROOT / "main"
SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
REF_DIR = REPO_ROOT / "reference"

# Define study-specific directories
STUDY_DIR = MAIN_DIR / study_name
RAW_DIR = STUDY_DIR / "raw_data"
QIIME_DIR = STUDY_DIR / "qiime2_artifacts"
PLOTS_DIR = STUDY_DIR / "plots"
TABLES_DIR = STUDY_DIR / "tables"
TMP_DIR = STUDY_DIR / 'tmp'

# In your main Snakefile, after defining directories:
log(f"[MAIN] SNAKEFILE_DIR = {SNAKEFILE_DIR}")
log(f"[MAIN] REPO_ROOT = {REPO_ROOT}")
log(f"[MAIN] STUDY_DIR = {STUDY_DIR}")

taxa_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Define metadata
metadata_path = STUDY_DIR / "metadata.tsv"
metadata_tsv = pd.read_csv(metadata_path, sep="\t", index_col=0)
log(f"[MAIN] Expected metadata path = {STUDY_DIR / 'metadata.tsv'}")
log(f"[MAIN] metadata.tsv exists? {os.path.exists(metadata_path)}")
log(f"[MAIN] Reading in .tsv from absolute path: {metadata_path.absolute()}...")

# The second column (first after the sample-ID index) is always the primary
# grouping axis, regardless of what it is named.
PRIMARY_GROUP_COL = metadata_tsv.columns[0]
log(f"[MAIN] Primary group column (2nd metadata column): '{PRIMARY_GROUP_COL}'")
if PRIMARY_GROUP_COL not in GROUPING_AXES:
    GROUPING_AXES.insert(0, PRIMARY_GROUP_COL)
    log(f"[MAIN] '{PRIMARY_GROUP_COL}' not in configured factors — auto-added as primary grouping axis.")



# ==============================
# Colorblind-friendly palette (Okabe-Ito)
# ==============================
OKABE_ITO = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
    "#999999", "#332288",
]

def _natural_sort_key(value):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", str(value))]

def _axis_order(df, col):
    """Return the display order for unique values of *col*.

    'PrimaryOrder' is reserved exclusively for the primary grouping column (the
    second column in the metadata file).  Secondary/tertiary axes always use
    natural sort — if order matters for those axes, move them to the second column.

    For the primary column:
      - If 'PrimaryOrder' is absent: natural sort.
      - If 'PrimaryOrder' is present and maps 1-to-1: use it.
      - If 'PrimaryOrder' is present but ambiguous: hard error — fix the metadata.
    """
    unique = df[col].dropna().unique()

    if col != PRIMARY_GROUP_COL or "PrimaryOrder" not in df.columns:
        return sorted(unique, key=_natural_sort_key)

    nunique_per_group = (
        df[[col, "PrimaryOrder"]]
        .dropna(subset=[col, "PrimaryOrder"])
        .groupby(col)["PrimaryOrder"]
        .nunique()
    )
    ambiguous = nunique_per_group[nunique_per_group > 1].index.tolist()
    if ambiguous:
        raise SystemExit(
            f"[ERROR] The 'PrimaryOrder' column does not map 1-to-1 with the primary grouping "
            f"column '{col}'. The following groups have conflicting PrimaryOrder values: {ambiguous}.\n"
            f"Please check your metadata — every row belonging to the same '{col}' value "
            f"must share the same PrimaryOrder integer.\n"
            f"Note: 'PrimaryOrder' is reserved for the primary grouping column (2nd metadata column). "
            f"If you intended to order a different axis, move that column to the second position."
        )

    return (
        df[[col, "PrimaryOrder"]]
        .dropna(subset=[col, "PrimaryOrder"])
        .groupby(col)["PrimaryOrder"]
        .first()
        .astype(int)
        .sort_values()
        .index.tolist()
    )

def _axis_colors(df, col, ordered_vals):
    """Return a {value: hex_color} mapping for *col*.

    'PrimaryColor' is reserved exclusively for the primary grouping column (the
    second column in the metadata file).  Secondary/tertiary axes always use
    Okabe-Ito — if custom colors matter for those axes, move them to the second
    column.

    For the primary column:
      - If 'PrimaryColor' is absent: Okabe-Ito.
      - If 'PrimaryColor' is present and maps 1-to-1 with full coverage: use it.
      - If 'PrimaryColor' is ambiguous or incomplete: hard error — fix the metadata.

    'PrimaryOrder' and 'PrimaryColor' are independent — either, both, or neither
    may be present without affecting the other.
    """
    if col != PRIMARY_GROUP_COL or "PrimaryColor" not in df.columns:
        return {val: OKABE_ITO[i % len(OKABE_ITO)] for i, val in enumerate(ordered_vals)}

    sub = df[[col, "PrimaryColor"]].dropna(subset=[col, "PrimaryColor"])

    nunique_per_group = sub.groupby(col)["PrimaryColor"].nunique()
    ambiguous = nunique_per_group[nunique_per_group > 1].index.tolist()
    if ambiguous:
        raise SystemExit(
            f"[ERROR] The 'PrimaryColor' column is ambiguous for the primary grouping "
            f"column '{col}'. The following groups have more than one color assigned: {ambiguous}.\n"
            f"Please check your metadata — every row belonging to the same '{col}' value "
            f"must share the same 'PrimaryColor' hex code.\n"
            f"Note: 'PrimaryColor' is reserved for the primary grouping column (2nd metadata column). "
            f"If you intended to color a different axis, move that column to the second position."
        )

    color_map = sub.groupby(col)["PrimaryColor"].first().to_dict()

    missing = [v for v in ordered_vals if v not in color_map]
    if missing:
        raise SystemExit(
            f"[ERROR] The 'PrimaryColor' column is declared but the following '{col}' groups "
            f"have no color entry: {missing}.\n"
            f"Please check your metadata — every value in the primary grouping column '{col}' "
            f"must have a 'PrimaryColor' hex code, or remove the 'PrimaryColor' column entirely "
            f"to fall back to the Okabe-Ito palette."
        )

    return {val: color_map[val] for val in ordered_vals}

GROUP_ORDERS = {axis: _axis_order(metadata_tsv, axis) for axis in GROUPING_AXES}
GROUP_COLORS = {
    axis: _axis_colors(metadata_tsv, axis, GROUP_ORDERS[axis])
    for axis in GROUPING_AXES
}

for _axis in GROUPING_AXES:
    if _axis != PRIMARY_GROUP_COL and len(GROUP_ORDERS[_axis]) > len(OKABE_ITO):
        log(f"[WARNING] Secondary axis '{_axis}' has {len(GROUP_ORDERS[_axis])} unique values but only {len(OKABE_ITO)} Okabe-Ito colors — colors will cycle.")

log(f"[MAIN] GROUP_ORDERS = {GROUP_ORDERS}")
log(f"[MAIN] GROUP_COLORS = {GROUP_COLORS}")



# Define second level of directories
TAXA_BARPLOT_DIR = PLOTS_DIR / "taxa_barplots"
DIVERSITY_DIR = PLOTS_DIR / "diversity"
ALPHA_DIR = DIVERSITY_DIR / "alpha_diversity"
BETA_DIR = DIVERSITY_DIR / "beta_diversity"
DIFFERENTIAL_ABUNDANCE_DIR = PLOTS_DIR / "differential_abundance"
CORE_METRICS_DIR = QIIME_DIR / "core-metrics-results"

# Define main qiime2 container
QIIME_CONTAINER = str(WORKFLOW_DIR / "containers/qiime2-2025.10-amplicon-core.sif")
# Named for the microbiomeMarker LEfSe that has since been removed; the image is
# kept because lefse_Huttenhower_run and every .res plot rule still need it, and
# renaming it would mean rebuilding the .sif for nothing. Its microbiomeMarker R
# package is now simply unused.
LEFSE_R_CONTAINER = str(WORKFLOW_DIR / "containers/lefse_MicrobiomeMarker-r-env.sif")
# The effect-size/cladogram renderer is Python and matplotlib and needs no R. This is
# simply the one image carrying the lefse package its two drawing scripts import.
EFFECTSIZE_CLADOGRAM_CONTAINER = LEFSE_R_CONTAINER

##############################################
# FINAL TARGETS
##############################################
# Define standard base qiime artifacts
dada2_outputs = [
    str(QIIME_DIR / "table-dada2.qzv"),
    str(QIIME_DIR / "rep-seqs-dada2.qzv"),
    str(QIIME_DIR / "dada2-stats.qzv"),
]

analysis_candidates_outputs = [
    str(QIIME_DIR / "table-analysis.qzv"),
    str(QIIME_DIR / "rep-seqs-analysis.qzv"),
    str(QIIME_DIR / "table-analysis.qzv"),
    str(QIIME_DIR / "rep-seqs-analysis.qzv")
]


taxa_barplot_outputs = [
    str(TAXA_BARPLOT_DIR / f"{db}" / f"taxa_barplot_{taxa_level}_by_{group}_{suffix}.svg")
    for db in reference_db
    for group in GROUPING_AXES
    for taxa_level in taxa_levels
    for suffix in ("samples", "groups")
]

taxa_barplot_png_outputs = [
    str(TAXA_BARPLOT_DIR / f"{db}" / f"taxa_barplot_{taxa_level}_by_{group}_{suffix}.png")
    for db in reference_db
    for group in GROUPING_AXES
    for taxa_level in taxa_levels
    for suffix in ("samples", "groups")
]


# Alpha diversity outputs, one standalone figure per metric (wildcards.metric).
# Shannon/Evenness/Simpson/Chao1 don't need a tree, so they're built for every
# db in reference_db; Faith PD does, so it's built only for
# phylogenetic_reference_db.
ALPHA_METRICS_NONPHYLO = ["shannon", "evenness", "simpson", "chao1"]
ALPHA_METRICS_PHYLO = ["faith_pd"]

alpha_metric_outputs = expand(
    ALPHA_DIR / "{db}_{metric}_{group}.svg",
    db=reference_db, metric=ALPHA_METRICS_NONPHYLO, group=GROUPING_AXES
) + expand(
    ALPHA_DIR / "{db}_{metric}_{group}.svg",
    db=phylogenetic_reference_db, metric=ALPHA_METRICS_PHYLO, group=GROUPING_AXES
)

alpha_metric_png_outputs = expand(
    ALPHA_DIR / "{db}_{metric}_{group}.png",
    db=reference_db, metric=ALPHA_METRICS_NONPHYLO, group=GROUPING_AXES
) + expand(
    ALPHA_DIR / "{db}_{metric}_{group}.png",
    db=phylogenetic_reference_db, metric=ALPHA_METRICS_PHYLO, group=GROUPING_AXES
)

alpha_sentinels = expand(
    DIVERSITY_DIR / ".{db}_{metric}_{group}_done",
    db=reference_db, metric=ALPHA_METRICS_NONPHYLO, group=GROUPING_AXES
) + expand(
    DIVERSITY_DIR / ".{db}_{metric}_{group}_done",
    db=phylogenetic_reference_db, metric=ALPHA_METRICS_PHYLO, group=GROUPING_AXES
)

# ------------------------------
# Build list of differential abundance outputs
# ------------------------------
differential_abundance_outputs = []
for db in reference_db:
    for method in DA_METHODS:
        
        # Check for any variation of LEfSe
        if "LEfSe" in method:
            # Determine which axes to use based on the method name
            if method == "LEfSe_per_factor":
                axes = DESIGN_INFO.get("factors", [])
            else: 
                axes = DESIGN_INFO.get("grouping_axes", [])

            for group in axes:
                if run_lefse_huttenhower:
                    differential_abundance_outputs.extend([
                        str(DIFFERENTIAL_ABUNDANCE_DIR / f"{db}/lefse_Huttenhower_LDA_by_{group}.svg"),
                        str(DIFFERENTIAL_ABUNDANCE_DIR / f"{db}/lefse_Huttenhower_Cladogram_by_{group}.svg"),
                    ])


differential_abundance_png_outputs = []
for db in reference_db:
    for method in DA_METHODS:
        if "LEfSe" in method:
            if method == "LEfSe_per_factor":
                axes = DESIGN_INFO.get("factors", [])
            else:
                axes = DESIGN_INFO.get("grouping_axes", [])
            for group in axes:
                if run_lefse_huttenhower:
                    differential_abundance_png_outputs.extend([
                        str(DIFFERENTIAL_ABUNDANCE_DIR / f"{db}/lefse_Huttenhower_LDA_by_{group}.png"),
                        str(DIFFERENTIAL_ABUNDANCE_DIR / f"{db}/lefse_Huttenhower_LDA_by_{group}.svg"),
                        str(DIFFERENTIAL_ABUNDANCE_DIR / f"{db}/lefse_Huttenhower_Cladogram_by_{group}.png"),
                        str(DIFFERENTIAL_ABUNDANCE_DIR / f"{db}/lefse_Huttenhower_Cladogram_by_{group}.svg"),
                    ])

picrust2_outputs = [
    str(STUDY_DIR / "picrust2_described" / "KO_metagenome_unstrat_described.tsv.gz"),
    str(STUDY_DIR / "picrust2_described" / "EC_metagenome_unstrat_described.tsv.gz"),
    str(STUDY_DIR / "picrust2_described" / "pathway_abun_unstrat_described.tsv.gz"),
    str(STUDY_DIR / "plots" / "picrust2_KO.png"),
    str(STUDY_DIR / "plots" / "picrust2_KO.svg"),
    str(STUDY_DIR / "plots" / "picrust2_EC.png"),
    str(STUDY_DIR / "plots" / "picrust2_EC.svg"),
    str(STUDY_DIR / "plots" / "picrust2_MetaCyc.png"),
    str(STUDY_DIR / "plots" / "picrust2_MetaCyc.svg"),
] if run_picrust2 else []

# Taxonomy artifacts (per DB)
taxonomy_outputs = expand(
    str(QIIME_DIR / "{db}-taxonomy.qza"),
    db=reference_db
) + expand(
    str(QIIME_DIR / "{db}-taxonomy.qzv"),
    db=reference_db
) + expand(
    str(QIIME_DIR / "{db}-taxa-bar-plots.qzv"),
    db=reference_db
)

# Non-phylogenetics alpha diversity vectors (computed once, DB-agnostic)
alpha_core_metrics_outputs = expand(
    str(CORE_METRICS_DIR / "{metric}-vector.qza"),
    metric=["shannon", "simpson", "evenness", "chao1"]
)

# Phylogenetic alpha diversity (per DB with a tree)
alpha_phylogenetic_outputs = expand(
    str(CORE_METRICS_DIR / "{db}-faith-pd-vector.qza"),
    db=phylogenetic_reference_db
)

# Non-phylogenetic PCoA (computed once, DB-agnostic)
pcoa_outputs = [
    str(CORE_METRICS_DIR / "jaccard-pcoa-results.qza"),
    str(CORE_METRICS_DIR / "bray-curtis-pcoa-results.qza"),
]

# Phylogenetic PCoA (per DB with a tree)
pcoa_phylogenetic_outputs = expand(
    str(CORE_METRICS_DIR / "{db}-unweighted-unifrac-pcoa-results.qza"),
    db=phylogenetic_reference_db
) + expand(
    str(CORE_METRICS_DIR / "{db}-weighted-unifrac-pcoa-results.qza"),
    db=phylogenetic_reference_db
)

# Beta diversity distance matrices (non-phylogenetic, computed once)
beta_distance_outputs = [
    str(CORE_METRICS_DIR / "jaccard-distance-matrix.qza"),
    str(CORE_METRICS_DIR / "bray-curtis-distance-matrix.qza"),
]

# Beta diversity distance matrices (phylogenetic, per DB with a tree)
beta_distance_phylogenetic_outputs = expand(
    str(CORE_METRICS_DIR / "{db}-unweighted-unifrac-distance-matrix.qza"),
    db=phylogenetic_reference_db
) + expand(
    str(CORE_METRICS_DIR / "{db}-weighted-unifrac-distance-matrix.qza"),
    db=phylogenetic_reference_db
)

# Beta diversity plots, one standalone figure per metric (wildcards.metric).
# Jaccard/Bray-Curtis don't need a tree, so they're built for every db in
# reference_db; the two UniFrac metrics do, so phylogenetic_reference_db only.
BETA_METRICS_NONPHYLO = ["jaccard", "braycurtis"]
BETA_METRICS_PHYLO = ["weighted_unifrac", "unweighted_unifrac"]

beta_metric_outputs = expand(
    str(BETA_DIR / "{db}_beta_{metric}_{group}.svg"),
    db=reference_db, metric=BETA_METRICS_NONPHYLO, group=GROUPING_AXES
) + expand(
    str(BETA_DIR / "{db}_beta_{metric}_{group}.svg"),
    db=phylogenetic_reference_db, metric=BETA_METRICS_PHYLO, group=GROUPING_AXES
)

beta_metric_png_outputs = expand(
    str(BETA_DIR / "{db}_beta_{metric}_{group}.png"),
    db=reference_db, metric=BETA_METRICS_NONPHYLO, group=GROUPING_AXES
) + expand(
    str(BETA_DIR / "{db}_beta_{metric}_{group}.png"),
    db=phylogenetic_reference_db, metric=BETA_METRICS_PHYLO, group=GROUPING_AXES
)

beta_sentinels = expand(
    str(DIVERSITY_DIR / ".{db}_beta_{metric}_{group}_done"),
    db=reference_db, metric=BETA_METRICS_NONPHYLO, group=GROUPING_AXES
) + expand(
    str(DIVERSITY_DIR / ".{db}_beta_{metric}_{group}_done"),
    db=phylogenetic_reference_db, metric=BETA_METRICS_PHYLO, group=GROUPING_AXES
)

# PERMANOVA (per DB — still runs for every db in reference_db, just over fewer
# distance matrices when that db has no tree; see permanova_dist_inputs)
permanova_outputs = expand(
    str(TABLES_DIR / "{db}_permanova_permdisp.tsv"),
    db=reference_db
)

# Report PDF: all visualisation plots collected into one document
report_plot_inputs = (
    taxa_barplot_outputs
    + list(alpha_metric_outputs)
    + list(beta_metric_outputs)
    + differential_abundance_outputs
    + picrust2_outputs
)

visualisations_pdf_output = str(STUDY_DIR / "visualisations_report.pdf")

print("taxa_barplot_outputs:", taxa_barplot_outputs)
print("dada2_outputs:", dada2_outputs)
print("analysis_candidates_outputs:", analysis_candidates_outputs)
print("taxonomy_outputs:", taxonomy_outputs)
print("alpha_core_metrics_outputs:", alpha_core_metrics_outputs)
print("beta_distance_outputs:", beta_distance_outputs)
print("pcoa_outputs:", pcoa_outputs)
print("alpha_metric_outputs:", alpha_metric_outputs)
print("beta_metric_outputs:", beta_metric_outputs)
print("permanova_outputs:", permanova_outputs)
print("picrust2_outputs:", picrust2_outputs)
print("differential_abundance_outputs:", differential_abundance_outputs)
print("alpha_phylogenetic_outputs:", alpha_phylogenetic_outputs)
print("beta_distance_phylogenetic_outputs:", beta_distance_phylogenetic_outputs)
print("pcoa_phylogenetic_outputs:", pcoa_phylogenetic_outputs)
print("alpha_sentinels:", alpha_sentinels)
print("beta_sentinels:", beta_sentinels)
print("visualisations_pdf_output:", visualisations_pdf_output)

rule all:
    input:
        *dada2_outputs,
        *analysis_candidates_outputs,
        *taxonomy_outputs,
        *taxa_barplot_outputs,
        *taxa_barplot_png_outputs,
        *alpha_core_metrics_outputs,
        *alpha_phylogenetic_outputs,
        *beta_distance_outputs,
        *beta_distance_phylogenetic_outputs,
        *pcoa_outputs,
        *pcoa_phylogenetic_outputs,
        *alpha_metric_outputs,
        *list(alpha_metric_png_outputs),
        *alpha_sentinels,
        *beta_metric_outputs,
        *list(beta_metric_png_outputs),
        *beta_sentinels,
        *permanova_outputs,
        *picrust2_outputs,
        *differential_abundance_outputs,
        *differential_abundance_png_outputs,
        visualisations_pdf_output,


rule pipeline_complete:
    input:
        *dada2_outputs,
        *analysis_candidates_outputs,
        *taxonomy_outputs,
        *taxa_barplot_outputs,
        *taxa_barplot_png_outputs,
        *alpha_core_metrics_outputs,
        *alpha_phylogenetic_outputs,
        *beta_distance_outputs,
        *beta_distance_phylogenetic_outputs,
        *pcoa_outputs,
        *pcoa_phylogenetic_outputs,
        *alpha_metric_outputs,
        *list(alpha_metric_png_outputs),
        *alpha_sentinels,
        *beta_metric_outputs,
        *list(beta_metric_png_outputs),
        *beta_sentinels,
        *permanova_outputs,
        *picrust2_outputs,
        *differential_abundance_outputs,
        *differential_abundance_png_outputs,
        visualisations_pdf_output,
    output:
        sentinel = str(STUDY_DIR / ".pipeline_complete")
    message:
        """
        [MAIN] Pipeline complete! All outputs have been generated successfully.
        """
    log:
        STUDY_DIR / "logs/pipeline_complete.log"
    shell:
        """
        exec > {log} 2>&1
        touch {output.sentinel}
        """

###############################################
# STEP 0 - Create directories
###############################################

rule make_dirs:
    output:
        marker = str(STUDY_DIR / ".dirs_created")
    run:
        STUDY_DIR.mkdir(parents=True, exist_ok=True)
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        TMP_DIR.mkdir(parents=True, exist_ok=True)
        QIIME_DIR.mkdir(parents=True, exist_ok=True)
        PLOTS_DIR.mkdir(parents=True, exist_ok=True)
        TABLES_DIR.mkdir(parents=True, exist_ok=True)
        TAXA_BARPLOT_DIR.mkdir(parents=True, exist_ok=True)
        DIVERSITY_DIR.mkdir(parents=True, exist_ok=True)
        ALPHA_DIR.mkdir(parents=True, exist_ok=True)
        BETA_DIR.mkdir(parents=True, exist_ok=True)
        DIFFERENTIAL_ABUNDANCE_DIR.mkdir(parents=True, exist_ok=True)
        CORE_METRICS_DIR.mkdir(parents=True, exist_ok=True)
        #  make dummy file to mark completion
        with open(output.marker, "w") as f:
            f.write("done")

###############################################
# STEP 1 — Build Manifest
###############################################

rule generate_manifest:
    input:
        make_dirs_marker = STUDY_DIR / ".dirs_created"
    output:
        manifest_file = STUDY_DIR / "manifest.tsv"
    params:
        raw_dir = RAW_DIR,
        sequence_type = sequence_type,
        region = region,
        study_dir = STUDY_DIR
    container:
        QIIME_CONTAINER
    log:
        STUDY_DIR / "logs/generate_manifest.log"
    script:
        "scripts/build_manifest.py"
        

###############################################
# STEP 2 - Trimming (Cutadapt) and Denoising (DADA2)
###############################################

include: f"rules/{sequence_type}.smk"

###############################################
# STEP 2.5 - Filter feature table to metadata samples
###############################################

rule filter_table_to_metadata:
    input:
        table    = QIIME_DIR / "table-dada2.qza",
        metadata = STUDY_DIR / "metadata.tsv"
    output:
        filtered = QIIME_DIR / "table-analysis.qza"
    container:
        QIIME_CONTAINER
    message:
        """[QIIME] Filtering feature table to analysis candidates...
        """
    log:
        STUDY_DIR / "logs/filter_table_to_metadata.log"
    shell:
        """
        exec > {log} 2>&1
        qiime feature-table filter-samples \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --o-filtered-table {output.filtered}
        """

rule visualise_table_for_analysis:
    input:
        table = QIIME_DIR / "table-analysis.qza"
    output:
        viz = QIIME_DIR / "table-analysis.qzv"
    container:
        QIIME_CONTAINER
    message:
        """
        [QIIME] Visualizing feature table of analysis candidates...
        """
    log:
        STUDY_DIR / "logs/visualise_table_for_analysis.log"
    shell:
        """
        exec > {log} 2>&1
        qiime feature-table summarize \
            --i-table {input.table} \
            --o-visualization {output.viz}
        """

rule filter_seqs_to_table:
    input:
        seqs  = QIIME_DIR / "rep-seqs-dada2.qza",
        table = QIIME_DIR / "table-analysis.qza"
    output:
        filtered = QIIME_DIR / "rep-seqs-analysis.qza"
    container:
        QIIME_CONTAINER
    message:
        """
        [QIIME] Filtering representative sequences to analysis candidates...
        """
    log:
        STUDY_DIR / "logs/filter_seqs_to_table.log"
    shell:
        """
        exec > {log} 2>&1
        qiime feature-table filter-seqs \
            --i-data {input.seqs} \
            --i-table {input.table} \
            --o-filtered-data {output.filtered}
        """

rule visualise_seqs_for_analysis:
    input:
        seqs = QIIME_DIR / "rep-seqs-analysis.qza"
    output:
        viz = QIIME_DIR / "rep-seqs-analysis.qzv"
    container:
        QIIME_CONTAINER
    message:
        """
        [QIIME] Visualizing representative sequences of analysis candidates...
        """
    log:
        STUDY_DIR / "logs/visualise_seqs_for_analysis.log"
    shell:
        """
        exec > {log} 2>&1
        qiime feature-table tabulate-seqs \
            --i-data {input.seqs} \
            --o-visualization {output.viz}
        """
###############################################
# STEP 3 - Taxonomy Classification and Phylogeny Construction
###############################################

for db in reference_db:
    include: f"rules/{db}.smk"

# Taxa barplot plotting is DB-agnostic (wildcarded on {db}) and must be included
# unconditionally — it cannot live inside a per-DB rules file like Greengenes2.smk
# or Silva138.smk, since only one of those may be included depending on REFERENCE_DB.
include: "rules/taxa_barplots.smk"

##############################################
# STEP 4 - Alpha and Beta Diversity Analyses
##############################################

include: "rules/diversity.smk"

##############################################
# STEP 5 - Differential Abundance Analyses
##############################################

DA_QIIME_TABLE = QIIME_DIR / "table-analysis.qza"
DA_QIIME_TAXA  = QIIME_DIR / "{db}-taxonomy.qza"

include: "rules/LEfSe.smk"

# ##############################################
# #STEP 6 - Functional Prediction with PICRUSt2
# ##############################################

include: "rules/PICRUSt2.smk"
include: "rules/run_audit.smk"

##############################################
# STEP 7 - Compile visualisations report PDF
##############################################

rule compile_visualisations_pdf:
    input:
        plots       = report_plot_inputs,
        metadata    = str(STUDY_DIR / "metadata.tsv"),
        permanova   = permanova_outputs,
        rarefy_base = str(TABLES_DIR / "rarefy_depth.csv"),
        rarefy_dbs  = expand(str(TABLES_DIR / "{db}_rarefy_depth.csv"), db=phylogenetic_reference_db),
    output:
        pdf = visualisations_pdf_output
    params:
        reference_dbs = reference_db,
        phylogenetic_reference_dbs = phylogenetic_reference_db,
        trunc_len_csv = str(TABLES_DIR / "trunc_len.csv") if sequence_type == "NGS" else None,
    container:
        str(WORKFLOW_DIR / "containers/report-env.sif")
    log:
        STUDY_DIR / "logs/compile_visualisations_pdf.log"
    script:
        "scripts/compile_visualisations_pdf.py"