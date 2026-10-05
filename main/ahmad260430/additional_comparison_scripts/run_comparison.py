#!/usr/bin/env python3
"""Orchestrate the Day17-vs-Day25 comparison, held within each treatment arm.

'Group' in metadata.tsv is the bare treatment name (e.g. 'IFA'), and the
separate 'Timepoint' column carries Day17/Day25 — so a straight
Day17-vs-Day25 comparison needs the treatment held fixed and grouping done on
'Timepoint'. Runs entirely as plain Python + subprocess calls into the
project's existing apptainer containers (qiime2 / LEfSe) — no Snakemake
involved.

Only two arms (IFA, IFA+IR) currently have both timepoints; Control,
Control+IR, Inosine and Inosine+IR only have Day17 samples so far and are
left out of TREATMENTS until Day25 data exists for them.

Only asks for outputs whose upstream artifacts already exist on disk for a
given reference db; it never kicks off a taxonomy/classification or
phylogeny/SEPP job itself. Run it again later (unchanged) once more artifacts
land and it will pick up whatever's newly available.

Usage (from anywhere):
    python3 run_comparison.py
    python3 run_comparison.py --only-treatment IFA
    python3 run_comparison.py --dry-run
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_subset_metadata import build_subset  # noqa: E402

THIS_DIR = Path(__file__).resolve().parent          # main/ahmad260430/day17_day25
STUDY_DIR = THIS_DIR.parent                          # main/ahmad260430
REPO_ROOT = STUDY_DIR.parent.parent                  # repo root

QIIME_DIR = STUDY_DIR / "qiime2_artifacts"
CORE_METRICS_DIR = QIIME_DIR / "core-metrics-results"
TABLES_DIR = STUDY_DIR / "tables"
METADATA = STUDY_DIR / "metadata.tsv"
LEFSE_PLUGIN_DIR = REPO_ROOT / "scripts" / "lefse"

QIIME_CONTAINER = REPO_ROOT / "containers" / "qiime2-2025.10-amplicon-core.sif"
LEFSE_CONTAINER = REPO_ROOT / "containers" / "lefse_MicrobiomeMarker-r-env.sif"
BIOM_R_CONTAINER = REPO_ROOT / "containers" / "ancombc2-env.sif"    # R with biomformat
MAASLIN3_CONTAINER = REPO_ROOT / "containers" / "maaslin3-env.sif"  # MaAsLin3 1.5.7

TREATMENTS = {
    "IFA":    "IFA",
    "IFA_IR": "IFA+IR",
}
FACTOR = "Timepoint"
DAY_ORDER = ["Day17", "Day25"]
DAY_COLORS = {"Day17": "#E69F00", "Day25": "#56B4E9"}  # Okabe-Ito, matches the main pipeline's palette

# Matches config/config_ahmad260430.yaml's LEfSe settings.
LEFSE_LDA_CUTOFF = "2.0"
LEFSE_KW_CUTOFF = "0.05"
LEFSE_WILCOXON_CUTOFF = "0.05"

COLLAPSE_DBS = ["Silva138", "Greengenes2"]
TAXA_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
ALPHA_METRICS_NONPHYLO = ["shannon", "evenness", "simpson", "chao1"]
BETA_METRICS_NONPHYLO = ["jaccard", "braycurtis"]
BETA_METRICS_PHYLO = ["weighted_unifrac", "unweighted_unifrac"]

DRY_RUN = False


OUT_DIR = STUDY_DIR / "comparison-day17-day25"

# One folder for the whole comparison, every artefact flat at its root under a
# self-describing name, `figures/` holding only the copies report.html embeds.
# Mirrors comparison-Treatment-x-Radiation-factorial, which is the other project-local
# analysis in this study; a reader who has seen one can navigate the other.
#
# The arm is a filename token rather than a directory because the two arms are
# one comparison asked twice, not two comparisons.
def art(kind, slug, db, *parts, ext):
    """<kind>_day17day25_<arm>[_<db>][_<part>...].<ext> at the comparison root.

    `db` is None for artefacts that do not read taxonomy. Shannon, Chao1,
    evenness, Simpson, Jaccard and Bray-Curtis come from database-independent
    QIIME artifacts, so naming a database in their filename would be a claim the
    file cannot support, and writing one copy per database would be the same
    figure twice.
    """
    tokens = [kind, "day17day25", slug] + ([db] if db else []) + [str(p) for p in parts]
    name = f"{'_'.join(tokens)}.{ext}"
    # Results in tables/, plots in figures/, prose at the root. Routing by what
    # the file is keeps the folder navigable however many parameterizations run.
    sub = OUT_DIR / ("figures" if ext in ("png", "svg")
                     else "tables" if ext in ("tsv", "res") else ".")
    sub.mkdir(parents=True, exist_ok=True)
    return sub / name


def color_args():
    args = []
    for day in DAY_ORDER:
        args += ["--color", f"{day}={DAY_COLORS[day]}"]
    return args


def run(cmd, label):
    print(f"[{label}] {' '.join(str(c) for c in cmd)}", flush=True)
    if DRY_RUN:
        return
    subprocess.run(cmd, check=True)


def apptainer_exec(container, cmd, label):
    full = ["apptainer", "exec", "--bind", f"{REPO_ROOT}:{REPO_ROOT}", str(container), *cmd]
    run(full, label)


def check_availability():
    def exists(p):
        return Path(p).exists()

    # Not per-database: the non-phylogenetic vectors and distance matrices carry
    # no {db} in their names and are the same file whichever database annotated
    # the study. Returning a list here is what produced one figure per database
    # from a single input.
    nonphylo_ready = exists(CORE_METRICS_DIR / "shannon-vector.qza")
    faith_pd_ready = [db for db in COLLAPSE_DBS if exists(CORE_METRICS_DIR / f"{db}-faith-pd-vector.qza")]
    unifrac_ready = [
        db for db in COLLAPSE_DBS
        if exists(CORE_METRICS_DIR / f"{db}-weighted-unifrac-pcoa-results.qza")
        and exists(CORE_METRICS_DIR / f"{db}-unweighted-unifrac-pcoa-results.qza")
    ]
    taxonomy_ready = [
        db for db in COLLAPSE_DBS
        if exists(QIIME_DIR / f"{db}-taxonomy.qza")
        and exists(TABLES_DIR / "exported-taxonomy" / f"{db}_taxonomy.tsv")
    ]

    for db in COLLAPSE_DBS:
        missing = []
        if db not in faith_pd_ready:
            missing.append("faith_pd")
        if db not in unifrac_ready:
            missing.append("UniFrac")
        if db not in taxonomy_ready:
            missing.append("LEfSe/taxa-barplots (no taxonomy)")
        if missing:
            print(f"[skip] {db}: {', '.join(missing)} — upstream artifact(s) not on disk yet.")

    return nonphylo_ready, faith_pd_ready, unifrac_ready, taxonomy_ready


def do_lefse(slug, db, subset_metadata, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    filtered_table = out_dir / f"table-analysis_{slug}.qza"

    apptainer_exec(QIIME_CONTAINER, [
        "qiime", "feature-table", "filter-samples",
        "--i-table", str(QIIME_DIR / "table-analysis.qza"),
        "--m-metadata-file", str(subset_metadata),
        "--o-filtered-table", str(filtered_table),
    ], f"{slug}/{db}/filter-table")

    lefse_res = art("LEfSe", slug, db, ext="res")

    apptainer_exec(LEFSE_CONTAINER, [
        "python3", str(THIS_DIR / "run_lefse_huttenhower.py"),
        "--feature-table", str(filtered_table),
        "--taxonomy", str(QIIME_DIR / f"{db}-taxonomy.qza"),
        "--metadata", str(subset_metadata),
        "--group-col", FACTOR,
        "--db", db,
        "--lda-cutoff", LEFSE_LDA_CUTOFF,
        "--kw-cutoff", LEFSE_KW_CUTOFF,
        "--wilcoxon-cutoff", LEFSE_WILCOXON_CUTOFF,
        *color_args(),
        "--lefse-dir", str(LEFSE_PLUGIN_DIR),
        "--output", str(lefse_res),
    ], f"{slug}/{db}/lefse-run")

    apptainer_exec(LEFSE_CONTAINER, [
        "python3", str(THIS_DIR / "plot_lefse_res.py"),
        "--res", str(lefse_res),
        "--db", db,
        # Same reason as the ANCOM-BC2 call: these land in the report as
        # arm-paired figures where left and right are named only in the caption.
        "--context", f"{TREATMENTS[slug]} arm",
        *color_args(),
        "--lda-svg", str(art("LEfSe", slug, db, "barplot", ext="svg")),
        "--cladogram-svg", str(art("LEfSe", slug, db, "cladogram", ext="svg")),
        "--lda-png", str(art("LEfSe", slug, db, "barplot", ext="png")),
        "--cladogram-png", str(art("LEfSe", slug, db, "cladogram", ext="png")),
    ], f"{slug}/{db}/lefse-plot")


def do_alpha(slug, db, metric, vector, subset_metadata, out_dir):
    apptainer_exec(QIIME_CONTAINER, [
        "python3", str(THIS_DIR / "plot_diversity.py"), "alpha",
        "--vector", str(vector),
        "--metadata", str(subset_metadata),
        "--metric", metric,
        "--factor", FACTOR,
        # Empty for a database-independent metric, which then prints no database
        # kicker rather than naming one it did not use.
        "--db", db or "",
        "--group-order", *DAY_ORDER,
        *color_args(),
        "--output-svg", str(art("alpha", slug, db, metric, ext="svg")),
        "--output-png", str(art("alpha", slug, db, metric, ext="png")),
    ], f"{slug}/{db or 'nodb'}/alpha-{metric}")


def do_beta(slug, db, metric, pcoa, subset_metadata, out_dir):
    apptainer_exec(QIIME_CONTAINER, [
        "python3", str(THIS_DIR / "plot_diversity.py"), "beta",
        "--pcoa", str(pcoa),
        "--metadata", str(subset_metadata),
        "--metric", metric,
        "--factor", FACTOR,
        "--db", db or "",
        "--group-order", *DAY_ORDER,
        *color_args(),
        "--output-svg", str(art("beta", slug, db, metric, ext="svg")),
        "--output-png", str(art("beta", slug, db, metric, ext="png")),
    ], f"{slug}/{db or 'nodb'}/beta-{metric}")


def do_taxa_barplot(slug, db, level, subset_metadata, out_dir):
    apptainer_exec(QIIME_CONTAINER, [
        "python3", str(THIS_DIR / "plot_taxa_barplot.py"),
        "--table-biom", str(TABLES_DIR / "study-seqs.biom"),
        "--taxonomy-tsv", str(TABLES_DIR / "exported-taxonomy" / f"{db}_taxonomy.tsv"),
        "--metadata-tsv", str(subset_metadata),
        "--db", db,
        "--level", level,
        "--factor", FACTOR,
        "--group-order", *DAY_ORDER,
        "--top-n", "20",
        "--arm", TREATMENTS[slug],
        # Rank over every arm's samples so the arms' figures share taxa and colours.
        "--rank-metadata-tsv", *(str(out_dir / "tables" / f"metadata_subset_{s}.tsv") for s in TREATMENTS),
        "--output-samples-svg", str(art("taxa", slug, db, level, "samples", ext="svg")),
        "--output-samples-png", str(art("taxa", slug, db, level, "samples", ext="png")),
        "--output-groups-svg", str(art("taxa", slug, db, level, "groups", ext="svg")),
        "--output-groups-png", str(art("taxa", slug, db, level, "groups", ext="png")),
    ], f"{slug}/{db}/taxa-{level}")


def do_report_inputs(db="Greengenes2"):
    """Everything report.html reads beyond the per-arm outputs, then a check of the page.

    Diversity statistics and their figures, the MaAsLin3 fits (both databases; the
    report uses Greengenes2), the structural-zero tables and heatmap, and the
    MaAsLin3 heatmaps and forest plots.
    """
    apptainer_exec(QIIME_CONTAINER, [
        "python3", str(THIS_DIR / "run_diversity_day17_day25.py"), str(OUT_DIR)], "diversity-stats")
    apptainer_exec(QIIME_CONTAINER, [
        "python3", str(THIS_DIR / "plot_beta_day17_day25.py")], "diversity-plots")
    for d in ("Greengenes2", "Silva138"):
        # ancombc2-env.sif only for its BIOM reader; nothing here fits ANCOM-BC2.
        apptainer_exec(BIOM_R_CONTAINER, [
            "Rscript", str(THIS_DIR / "export_counts_day17_day25.R"), "--db", d], f"{d}/export-counts")
        apptainer_exec(MAASLIN3_CONTAINER, [
            "Rscript", str(THIS_DIR / "run_maaslin3_day17_day25.R"), "--db", d], f"{d}/maaslin3")
    for script, args in (("structural_zeros_day17_day25.py", []),
                         ("plot_presence_day17_day25.py", []),
                         ("plot_maaslin3_day17_day25.py", []),
                         ("check_report_day17_day25.py", [])):
        apptainer_exec(QIIME_CONTAINER, ["python3", str(THIS_DIR / script), *args], script)


def main():
    global DRY_RUN
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only-treatment", choices=sorted(TREATMENTS))
    parser.add_argument("--dry-run", action="store_true", help="print the commands without running them")
    args = parser.parse_args()
    DRY_RUN = args.dry_run

    nonphylo_ready, faith_pd_ready, unifrac_ready, taxonomy_ready = check_availability()

    slugs = [args.only_treatment] if args.only_treatment else list(TREATMENTS)
    for slug in slugs:
        treatment_name = TREATMENTS[slug]
        out_dir = OUT_DIR
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "tables").mkdir(parents=True, exist_ok=True)
        subset_metadata = out_dir / "tables" / f"metadata_subset_{slug}.tsv"

        print(f"\n=== {slug} ({treatment_name!r}) ===")
        build_subset(METADATA, treatment_name, subset_metadata)

        for db in taxonomy_ready:
            do_lefse(slug, db, subset_metadata, out_dir)
            for level in TAXA_LEVELS:
                do_taxa_barplot(slug, db, level, subset_metadata, out_dir)

        # Database-independent, so they run once with db = None rather than once
        # per database. The artifacts these read carry no {db} in their names,
        # and taxonomy does not enter a Shannon index or a Bray-Curtis distance.
        if nonphylo_ready:
            for metric in ALPHA_METRICS_NONPHYLO:
                do_alpha(slug, None, metric, CORE_METRICS_DIR / f"{metric}-vector.qza", subset_metadata, out_dir)
            for metric in BETA_METRICS_NONPHYLO:
                pcoa_name = "jaccard-pcoa-results.qza" if metric == "jaccard" else "bray-curtis-pcoa-results.qza"
                do_beta(slug, None, metric, CORE_METRICS_DIR / pcoa_name, subset_metadata, out_dir)

        for db in faith_pd_ready:
            do_alpha(slug, db, "faith_pd", CORE_METRICS_DIR / f"{db}-faith-pd-vector.qza", subset_metadata, out_dir)

        for db in unifrac_ready:
            for metric in BETA_METRICS_PHYLO:
                pcoa_name = f"{db}-{metric.replace('_', '-')}-pcoa-results.qza"
                do_beta(slug, db, metric, CORE_METRICS_DIR / pcoa_name, subset_metadata, out_dir)

    # Both arms are fitted together, so this runs once, after every arm.
    if not args.only_treatment and "Greengenes2" in taxonomy_ready:
        do_report_inputs()

    print("\nDone.")


if __name__ == "__main__":
    main()
