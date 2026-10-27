import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

# ── Snakemake bindings ────────────────────────────────────────────────────────
feature_table  = snakemake.input.feature_table
taxonomy       = snakemake.input.taxonomy
metadata       = snakemake.input.metadata

lefse_results  = snakemake.output.lefse_Huttenhower_results

group_col       = snakemake.params.group_col
lda_cutoff      = snakemake.params.lda_cutoff
kw_cutoff       = snakemake.params.kw_cutoff
wilcoxon_cutoff = snakemake.params.wilcoxon_cutoff
colors_map      = snakemake.params.colors
db              = snakemake.wildcards.db
region          = snakemake.params.region

# ── Per-database taxonomic resolution ──────────────────────────────────────
# Greengenes2's full-length 16S classifier resolves real species-level calls
# (e.g. "s__Duncaniella muris"); Silva138's species field is empty for every
# ASV, so collapsing past genus would only add noise. A V3V4 amplicon is too
# short for a species call to be trusted, so Greengenes2 stops at genus there
# too. `collapse_level` is the QIIME 2 `taxa collapse` rank (6 = genus,
# 7 = species) used when building the LEfSe input table.
COLLAPSE_LEVEL = {
    "Silva138":    6,   # genus
    "Greengenes2": 7 if region == "full_length" else 6,   # species on full length only
}
if db not in COLLAPSE_LEVEL:
    raise ValueError(f"No LEfSe collapse level for db={db!r}. Known: {sorted(COLLAPSE_LEVEL)}")
collapse_level = COLLAPSE_LEVEL[db]

# ── Resolve the conda env's own Python ───────────────────────────────────────
# sys.executable is Snakemake's Python; we need the one inside the rule's
# conda env, which is always at $CONDA_PREFIX/bin/python.
CONDA_PYTHON = Path(sys.executable).parent / "python"
if not CONDA_PYTHON.exists():
    raise FileNotFoundError(
        f"Could not find conda env Python at {CONDA_PYTHON}.\n"
        f"sys.executable was: {sys.executable}"
    )

# ── Resolve lefse/ folder relative to this script ────────────────────────────
LEFSE_DIR = Path(snakemake.scriptdir) / "lefse"

def lefse(name: str) -> str:
    p = LEFSE_DIR / name
    if not p.exists():
        raise FileNotFoundError(f"LEfSe script not found: {p}")
    return str(p)

# ── Helper ────────────────────────────────────────────────────────────────────
def run(cmd: list) -> None:
    print("Running:", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True)

# ── A transient sort-key prefix on the class labels ───────────────────────────
# LEfSe assigns its internal class ordering from the class labels themselves.
# The figures need that ordering to match the configured group order, and the
# renderer reproduces this same prefix on a temp copy of the .res at draw time.
# Applying it here as well keeps LEfSe's own class ordering identical to what it
# has always been, so the LDA values this file records do not shift.
order = list(colors_map.keys())

# Both the original (space) and QIIME 2 (underscore) forms of each group name
# map to an underscore-only prefixed label. LEfSe scripts split class names on
# whitespace internally, so spaces in class labels cause parse errors;
# underscores are safe.
sortkey = {}
for i, g in enumerate(order):
    g_norm = g.replace(" ", "_")
    label  = f"zzsort{i:02d}_{g_norm}"
    sortkey[g]      = label
    sortkey[g_norm] = label

def strip_prefix(path) -> None:
    """Clean up the persisted .res file (a real pipeline output, unlike the
    SVG/PNG figures) so anything reading it downstream sees plain group names."""
    p = Path(path)
    text = p.read_text()
    for g in order:
        text = text.replace(sortkey[g], g)
    p.write_text(text)

# ── Intermediates as tempfiles ────────────────────────────────────────────────
with tempfile.TemporaryDirectory() as _tmp:
    tmpdir = Path(_tmp)
    tmpdir.mkdir(parents=True, exist_ok=True)
    lefse_input_tsv = tmpdir / "lefse_input.tsv"
    lefse_input_in  = tmpdir / "lefse_input.in"

    # Step 1 – QIIME 2 artefacts → LEfSe TSV
    run([
        CONDA_PYTHON, lefse("plugin_lefse_input.py"),
        "--table_file",    feature_table,
        "--taxonomy_file", taxonomy,
        "--metadata_file", metadata,
        "--output_file",   lefse_input_tsv,
        "--class_col",     group_col,
        "--level",         str(collapse_level),
    ])

    # Step 1.5 – rename class labels in the TSV group row with the sort-key
    # prefix so alphabetical order == desired order. No column reordering needed:
    # the plotting scripts re-derive class order from the labels themselves.
    tsv = pd.read_csv(lefse_input_tsv, sep="\t", header=None,
                      dtype=str, keep_default_na=False)
    feature_col = tsv.iloc[:, 0]          # feature-name column (untouched)
    data = tsv.iloc[:, 1:].copy()         # sample columns; row 0 = group labels

    present = set(data.iloc[0].unique())
    missing = present - set(sortkey)
    if missing:
        raise ValueError(
            f"Group(s) in the data are absent from colors_map: {sorted(missing)}.\n"
            f"colors_map keys: {sorted(sortkey)}"
        )

    data.iloc[0] = data.iloc[0].map(sortkey)
    pd.concat([feature_col, data], axis=1).to_csv(
        lefse_input_tsv, sep="\t", header=False, index=False
    )

    # Step 2 – format TSV → LEfSe .in
    run([
        CONDA_PYTHON, lefse("plugin_lefse_format.py"),
        lefse_input_tsv,
        lefse_input_in,
        "-c", "1",
        "-o", "1000000",
    ])

    # Step 3 – run LEfSe → results (class column now carries the prefixed labels)
    run([
        CONDA_PYTHON, lefse("plugin_lefse_run.py"),
        lefse_input_in,
        lefse_results,
        "-l", str(lda_cutoff),
        "-a", str(kw_cutoff),
        "-w", str(wilcoxon_cutoff),
    ])

    # Step 5 – restore original group names (spaces included) in the persisted
    # .res file so anything downstream sees clean labels.
    strip_prefix(lefse_results)
