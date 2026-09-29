# 16S rRNA Amplicon Sequencing Pipeline — Overview

This Snakemake pipeline processes 16S rRNA amplicon sequencing data (NGS or TGS) from raw reads through to diversity analyses, differential abundance testing, and functional prediction. Each stage is described below.

> **New to this repo?** See [`quickstart.md`](quickstart.md) for `metadata.tsv` column-naming rules, `config_project.yaml` setup, and how to actually invoke Snakemake (dry runs, `--until`, `--touch`, etc). This document assumes that setup is done and walks through what each stage of the workflow produces.

> **Plot outputs:** Every plot the pipeline generates is written in **two formats — one `.svg` and one `.png`** (e.g. `*.svg` and `*.png` of the same name).

---

## Stage 0 — Setup & Validation

**Rule:** `make_dirs`

Before any analysis begins, the pipeline:

- Validates the configuration (sequence type, targeted region, reference databases, experimental design).
- Creates all required output directories (`raw_data/`, `qiime2_artifacts/`, `plots/`, `tables/`, `tmp/`, and subdirectories for diversity, taxa barplots, and differential abundance).
- Reads the sample metadata (`metadata.tsv`) and determines which grouping axes (experimental factors) will drive downstream comparisons.

---

## Stage 1 — Manifest Building

**Rule:** `generate_manifest` → `scripts/build_manifest.py`

Scans the `raw_data/` directory and builds a QIIME2-compatible manifest TSV file that maps sample IDs to their raw FASTQ file paths. The manifest format differs slightly between NGS (paired-end) and TGS (single-end CCS) input.
Currently we work with Illumina for our NGS data and PacBio for our TGS data.

---

## Stage 2 — Adapter Trimming & Denoising

**Rules:** `NGS.smk` or `TGS.smk` (selected based on `SEQUENCE_TYPE`)

### 2a. Import

Raw reads are imported into QIIME2 as a demultiplexed artifact (`demux.qza`).

- **NGS**: imported as paired-end sequences (`PairedEndSequencesWithQuality`).
- **TGS**: imported as single-end sequences (`SequencesWithQuality`).

### 2b. Primer Trimming (Cutadapt)

- **NGS**: Adapter/primer sequences are removed using `qiime cutadapt trim-paired` with the standard V3V4 primer pair (341F/805R with Nextera adapter overhangs) at up to 10% error rate.
- **TGS**: The cutadapt step is a no-op — it renames `demux.qza` → `trimmed-demux.qza` without invoking cutadapt (CCS reads are already primer-free).

### 2c. Quality Profiling & Truncation Length Selection (NGS only)

- Quality scores are exported from the trimmed demux visualisation.
- `generate_trunc_length.py` automatically determines forward and reverse truncation lengths by finding where median quality drops below Q20 or when more than 10% of reads are being dropped.

### 2d. DADA2 Denoising

Sequences are denoised to produce **Amplicon Sequence Variants (ASVs)**:

- **NGS** (`denoise-paired`): paired-end reads are quality-filtered, denoised, merged, and chimera-filtered. Truncation lengths are applied automatically from Stage 2c.
- **TGS** (`denoise-ccs`): PacBio CCS reads are processed using the full-length denoising workflow with a 5′ primer anchor (`--p-front AGAGTTTGATCMTGGCTCAG`, the universal 27F primer). Truncation length is hard-coded to `0` (no truncation), which is correct for full-length CCS reads.

Outputs:

- `table-dada2.qza` — ASV feature table (samples × ASVs).
- `rep-seqs-dada2.qza` — representative sequences for each ASV.
- `dada2-stats.qza/.qzv` — per-sample denoising statistics.
- `base-transition-stats-dada2.qza` — base substitution/transition statistics from DADA2 (both NGS and TGS).

### 2e. Sample Filtering

After denoising, the feature table and representative sequences are automatically filtered to retain only samples that are present in `metadata.tsv`. This is the point at which you can edit `metadata.tsv` to exclude specific samples from all downstream analyses (e.g. failed libraries or QC outliers) — samples absent from the metadata are silently dropped. The result is:

- `table-analysis.qza` — filtered ASV feature table.
- `rep-seqs-analysis.qza` — filtered representative sequences.

### 2f. Rarefaction Depth Determination (Base)

The filtered feature table (`table-analysis.qza`) is exported and a base rarefaction depth is computed from per-sample read counts (`rarefy_depth.csv`). A configurable percentile (`RAREFY_DEPTH_PERCENTILE`) controls the selection: `0` selects the **minimum** read count across all samples; any other value (e.g. `10`) selects the corresponding percentile, which can preserve more samples when one sample has an unusually low count.

---

## Stage 3 — Taxonomy Classification & Phylogeny

**Rules:** `Greengenes2.smk` and/or `Silva138.smk` (one or both, based on `REFERENCE_DB`)

### 3a. Taxonomic Classification

Representative sequences are classified against a pre-trained Naive Bayes classifier for the chosen reference database and amplicon region:

| Database    | V3V4 classifier                               | Full-length classifier                                                 |
| ----------- | --------------------------------------------- | ---------------------------------------------------------------------- |
| Greengenes2 | `Greengenes2/2024.09.custom.V3V4.nb.qza`    | `Greengenes2/2024.09.backbone.full-length.nb.qza`                    |
| SILVA 138   | `Silva138/silva138_AB_V3-V4_classifier.qza` | `Silva138/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza` |

Classifier files live in `REF_DIR/{database}/`. This produces a taxonomy artifact (`.qza`) and a tabular visualisation (`.qzv`).

**SILVA 138 V3V4 runs in a pinned legacy container.** `silva138_AB_V3-V4_classifier.qza` was trained under QIIME 2 2021.4 (scikit-learn 0.24.1) with extract-reads primers matching `NGS.smk`'s cutadapt primers exactly. `classify-sklearn` in the current container (scikit-learn 1.4.2) refuses to unpickle a model trained on a different scikit-learn release rather than risk silent corruption, and no official SILVA 138 classifier matches this primer set at 1.4.2. Retraining would discard an already-validated model, so `rules/Silva138.smk` pins `containers/qiime2-2021.4-silva138-classify.sif` for this one rule instead.

That container also rejects artifacts written by the modern one, so `Silva138_reimport_rep_seqs_legacy` first re-wraps `rep-seqs-analysis.qza` into a 2021.4-native artifact (`rep-seqs-analysis-legacy.qza`) by pulling `dna-sequences.fasta` out of the zip with stdlib `zipfile` and re-importing it. Only the reads need this round trip; the classifier already loads natively. Greengenes2 and SILVA 138 full-length classify in the normal container.

### 3b. Genus-Level Collapse (V3V4 NGS only)

For V3V4 NGS data, the feature table is additionally collapsed to genus level (rank 6) using `qiime taxa collapse`. This produces `{db}-genus-table-collapsed.qza` and a corresponding exported `.biom` file, one per database (`{db}-study-seqs-genus.biom`). Differential abundance analyses in Stage 5 receive this genus-collapsed table rather than the ASV-level table, as genus is the finest reliable taxonomic resolution for the V3V4 region. Full-length TGS data passes the ASV-level table to Stage 5 instead — but note that two of the three DA methods then collapse it themselves, per database, before testing (see Stage 5).

### 3c. Phylogeny — Backbone Mapping

ASVs are mapped onto the full-length backbone sequence database using `qiime greengenes2 non-v4-16s` (or the equivalent for SILVA138). This produces:

- A **phylogeny-aware feature table** (`{db}-table.qza`) — ASVs remapped to backbone representatives.
- **Phylogeny-aware representative sequences** (`{db}-rep-seqs.qza`).

These are the inputs for phylogenetic diversity analyses. The actual **phylogenetic tree** used for UniFrac metrics (weighted and unweighted) is a **pre-built backbone `.nwk.qza`** that must already exist in `REF_DIR/{db}/`. The pipeline does not construct this tree itself — **except for Silva138**, which ships no pre-built backbone tree at all. For Silva138 the tree is built per-study via fragment-insertion (`qiime fragment-insertion sepp`, `rules/Silva138.smk`), placing ASVs onto QIIME2's SILVA 128 reference (no official SILVA 138 SEPP reference exists; taxonomy assignment is still the real SILVA 138 release). SEPP is memory-heavy and has previously OOM-killed the host, so:

- The rule logs a compact one-line-per-minute progress summary (subset counts, resource usage) via `scripts/sepp_monitor_lib.sh`; a job already in flight can also be attached to from any terminal with `scripts/watch_sepp.sh`.
- Setting `SILVA138_SKIP_PHYLOGENY: true` in the config skips building this tree entirely. When set, Silva138 is dropped from `phylogenetic_reference_db` (a narrower list than `reference_db`) and every target that needs a tree — Faith's PD, UniFrac, their PCoA, and the per-DB rarefaction depth — is built without it; taxonomy, taxa barplots, and non-phylogenetic diversity metrics are unaffected.

### 3d. Taxa Barplot Visualisations

- **QIIME interactive barplot** (`.qzv`) — generated natively in QIIME2 from the filtered ASV table.
- **Custom static barplots** (`.svg` + `.png`) — generated per taxonomic level (Kingdom → Species) and per grouping axis using `plot_taxa_barplot_levels.py`, showing the top 20 taxa per group.

For every taxonomic level × grouping axis combination, **two versions are always produced**, distinguished by the `_samples` vs `_groups` suffix — for example:

- `taxa_barplot_Class_by_Group_samples.png` — **per-sample** view: each sample gets its own bin (bar).
- `taxa_barplot_Class_by_Group_groups.png` — **per-group** view: each group is collapsed into a single bin based on pooled abundance across its samples.

(Each version is written as both `.svg` and `.png`, following the pipeline-wide plot output convention.)

### 3e. Artifact Export

The feature table (`.biom`), representative sequences (`.fna`), and taxonomy assignments (`.tsv`) are exported from QIIME2 for use in downstream non-QIIME tools (differential abundance, PICRUSt2).

---

## Stage 4 — Diversity Analyses

**Rules:** `NGS.smk` / `TGS.smk` (rarefaction), `diversity.smk`

### 4a. Two-Stage Rarefaction

The pipeline performs **two independent rarefaction steps**, each targeting a different feature table:

**Base rarefaction** (ASV-level, non-phylogenetic):
The filtered ASV table (`table-analysis.qza`) is rarefied to the depth computed in Stage 2f (`rarefy_depth.csv`), producing `table-analysis-rarefied.qza`. This rarefied table is used for all non-phylogenetic diversity metrics (Shannon, Simpson, Pielou's evenness, Jaccard, Bray-Curtis).

**Per-DB rarefaction** (backbone-mapped, phylogenetic):
Each reference database's backbone-mapped table (`{db}-table.qza`) is rarefied independently. A separate rarefaction depth is calculated per database by reading the backbone table's per-sample counts and applying the same percentile logic (`{db}_rarefy_depth.csv`). This produces `{db}-table-rarefied.qza` and is used for all phylogenetic diversity metrics (Faith's PD, weighted and unweighted UniFrac). The two depths may differ because backbone mapping can change per-sample read counts.

### 4b. Alpha Diversity

Five within-sample diversity metrics are calculated for each reference database — four non-phylogenetic, plus Faith's PD for databases in `phylogenetic_reference_db` (see Silva138's SEPP exception in 3c):

| Metric            | What it measures                           | Phylogenetic? |
| ----------------- | ------------------------------------------ | ------------- |
| Shannon index     | Species richness and evenness              | No            |
| Simpson index     | Dominance / evenness                       | No            |
| Pielou's evenness | How evenly abundances are distributed      | No            |
| Chao1             | Estimated richness, weighting rare taxa    | No            |
| Faith's PD        | Phylogenetic diversity (branch length sum) | Yes (per DB)  |

Each metric is plotted as its own standalone grouped box/strip plot per grouping axis (`plot_alpha_diversity_metric.py`, one figure per `{db, metric, group}` combination), written as both `.svg` and `.png`.

### 4c. Beta Diversity

Four between-sample dissimilarity matrices are computed (the two UniFrac metrics only for databases in `phylogenetic_reference_db`):

| Metric             | Phylogenetic? | What it captures                         | Input table     |
| ------------------ | ------------- | ---------------------------------------- | --------------- |
| Jaccard            | No            | Presence/absence overlap                 | Base rarefied   |
| Bray-Curtis        | No            | Abundance-weighted overlap               | Base rarefied   |
| Unweighted UniFrac | Yes           | Presence/absence + phylogenetic distance | Per-DB rarefied |
| Weighted UniFrac   | Yes           | Abundance-weighted phylogenetic distance | Per-DB rarefied |

PCoA ordinations are computed for all four matrices; each metric is plotted as its own standalone 2D scatter plot coloured by grouping axis (`plot_beta_diversity_metric.py`, one figure per `{db, metric, group}` combination), written as both `.svg` and `.png`.

### 4d. Statistical Testing (PERMANOVA / PermDisp)

For all four distance matrices and all grouping axes, the pipeline runs:

- **PERMANOVA** — tests whether group centroids differ significantly.
- **PermDisp** — tests whether group dispersions (spread) are homogeneous (a prerequisite for interpreting PERMANOVA).

Results are saved to `{db}_permanova_permdisp.tsv`.

---

## Stage 5 — Differential Abundance Analysis

**Rules:** `LEfSe.smk`, `ANCOMBC2.smk`

### Taxonomic resolution of the input table

Both methods start from the exported `.biom` feature table and taxonomy TSV. Which table that is depends on the region: V3V4 NGS gets the genus-collapsed table from Stage 3b, full-length TGS gets the ASV-level table.

Both methods then collapse further, **by database rather than by sequence type**, because resolution is a property of the classifier: Greengenes2's full-length classifier resolves real species, while Silva138's species field is empty for every ASV, so collapsing past genus would only add noise. Collapsing also makes the methods comparable in the shared `.res` renderer, whose nodes are taxonomy strings — uncollapsed ASVs sharing a taxonomy merge onto one node no matter how many were significant.

| Method                | Collapse applied inside the method                          |
| --------------------- | ----------------------------------------------------------- |
| `lefse_Huttenhower` | Silva138 → rank 6 (genus); Greengenes2 → rank 7 (species) |
| `ancombc2`          | same                                                        |

So on full-length TGS data both operate at genus/species level, not ASV level.

### LEfSe (Linear Discriminant Analysis Effect Size)

LEfSe is run with the original Huttenhower implementation, switchable via config (`RUN_LEFSE_HUTTENHOWER`) and producing its outputs as both `.svg` and `.png` when enabled.

Each method is split into a `_run` rule and a `_plot` rule. The `_run` rule executes that method's own statistics in that method's own environment and writes a **`.res` file**; the `_plot` rule renders it. The `.res` file is the only thing that crosses between them, so re-plotting never re-runs the statistics.

A `.res` file is a prefix-closed table over the taxonomy tree, one TAB-separated row per node: dot-separated lineage, `log10` of the largest per-class mean abundance, and then the enrichment class, LDA score and p-value for nodes the method called significant. Because every method emits this one shape, they are all rendered by the **same two drawing scripts** (`scripts/plotters/effectsize_barplot.py` and `scripts/plotters/cladogram.py`) via the shared `scripts/effectsize_cladogram_render.py`, so the figures are directly comparable across methods rather than each library drawing in its own idiom. Those scripts and the renderer are named for the plots they produce, not for LEfSe: the format originated there, but ANCOM-BC2 renders through the same path. Group colors and legend order come from `GROUP_COLORS` and are pinned to the full configured group universe, so a group keeps its color even in comparisons it does not appear in.

When a method finds nothing significant, the upstream behaviour is preserved: the LDA bar chart is written as an empty file and the cladogram renders the uncoloured background tree.

**`lefse_Huttenhower` — vanilla Python-style Huttenhower LEfSe** (`run_lefse_Huttenhower.py`, toggled by `RUN_LEFSE_HUTTENHOWER`, default `true`) — the original LEfSe implementation, which tests every taxonomic rank. Its results file is in `.res` format but carries a `.txt` extension (`lefse_Huttenhower_results_by_{group}.txt`), unlike the other two methods, which use `.res`. For each grouping axis it outputs:

- **`lefse_Huttenhower_LDA_by_{group}`** — LDA bar chart of ranked differentially abundant taxa.
- **`lefse_Huttenhower_Cladogram_by_{group}`** — cladogram of the enriched taxa.

### ANCOMBC2

Fits one compositionality-corrected model per grouping axis (`fix_formula` = the grouping column plus any configured covariates/interactions from `DESIGN`), toggled by `RUN_ANCOMBC2`. Which group a taxon is enriched in, and its effect size, come from the strongest pairwise contrast that clears mdFDR correction (Holm, bootstrap) — this correction is what distinguishes it from LEfSe's uncorrected pairwise tests.

How a taxon qualifies depends on how many levels the grouping axis has:

- **Three or more levels:** the global omnibus test runs, and a taxon must clear its FDR-adjusted q-value before any contrast is considered. Omnibus-significant taxa whose pairwise family localises nothing are dropped, since there is no single group to attribute them to.
- **Two levels:** `ancombc2()` deactivates the multi-group comparisons and returns no omnibus table, so there is no gate to apply. Every tested taxon is considered and the pairwise q-values are the only FDR control in play. 

Because a metadata column name cannot contain a space in ANCOM-BC2's `fix_formula` (it resolves variable names by parsing that string, and does not strip backticks), `run_ancombc2.R` builds a spaceless alias column for every name entering the model and reads the term names off the alias. Group *values* are unaffected and keep their spaces throughout.

`ancombc2_run` writes a results `.tsv` and a `.res`. The `.tsv` is the raw `ancombc2()` table, kept for inspection; only the `.res` is plotted.

- **`ancombc2_effectsize_cladogram_plot`** — the FDR-controlled hits rendered through the shared `.res` renderer (`ANCOMBC2_EffectSize_by_{group}`, `ANCOMBC2_Cladogram_by_{group}`), so ANCOM-BC2's stricter hit list can be visually compared against LEfSe's more permissive one in the same chart format. The bars are `|log2FC| (mdFDR-controlled)`, not LDA scores, which is why neither the rule nor the file says LDA. `ancombc2()` returns a natural-log fold change, so `run_ancombc2.R` divides by `log(2)` before writing the `.res`. Bar lengths are therefore not comparable with the LEfSe charts, whose scores are `log10(1 + effect)` over CPM-normalised abundances; only the taxa and their direction are.

---

## Stage 6 — Functional Prediction with PICRUSt2

**Rules:** `PICRUSt2.smk`

PICRUSt2 is toggled by `RUN_PICRUST2` (default `true` if unset — set explicitly to `false` in this project's config to skip it).

PICRUSt2 predicts the functional potential of the microbial community from 16S marker gene data alone, without requiring shotgun metagenomics.

### 6a. PICRUSt2 Pipeline

Using the exported `.fna` representative sequences and `.biom` feature table, the full PICRUSt2 pipeline is run to predict:

- **KEGG Orthology (KO)** metagenome abundances.
- **Enzyme Commission (EC)** number abundances.
- **MetaCyc pathway** abundances.

### 6b. Add Descriptions

PICRUSt2's own `add_descriptions.py` CLI tool (bundled in the PICRUSt2 environment, not a script in `scripts/`) annotates each KO, EC, and pathway ID with a human-readable description. Rule: `picrust2_describe`, writing to `picrust2_described/`.

### 6c. Visualisation

`plot_picrust2.R` generates three separate heatmaps of predicted abundances grouped by sample metadata — one each for KEGG Orthology, Enzyme Commission, and MetaCyc pathways (`picrust2_KO`, `picrust2_EC`, `picrust2_MetaCyc`), each written as both `.svg` and `.png`.

---

## Stage 7 — Visualisations Report

**Rule:** `compile_visualisations_pdf`

After all plots are generated, the pipeline compiles every visualisation into a single PDF report (`visualisations_report.pdf`) in the study directory. It takes whatever Stage 3-6 actually produced: taxa barplots, the per-metric alpha and beta diversity plots, every enabled differential abundance method (LEfSe **and** ANCOM-BC2), and the PICRUSt2 heatmaps. PERMANOVA/PermDisp results and rarefaction depth information are embedded as summary tables.

---

## Output Summary

| Stage                | Key outputs                                                                                                                                                                                                                                                                                   |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0 — Setup           | Directory structure, validated config                                                                                                                                                                                                                                                         |
| 1 — Manifest        | `manifest.tsv`                                                                                                                                                                                                                                                                              |
| 2 — Denoising       | `table-dada2.qza`, `rep-seqs-dada2.qza`, `dada2-stats.qzv`, `table-analysis.qza`, `rep-seqs-analysis.qza`                                                                                                                                                                           |
| 3 — Taxonomy        | `{db}-taxonomy.qza`, `{db}-taxa-bar-plots.qzv`, `taxa_barplot_*_{samples,groups}.{svg,png}`; V3V4 NGS also `{db}-study-seqs-genus.biom`                                                                                                                                               |
| 4 — Diversity       | Alpha/beta plots (`.svg` + `.png`), `{db}_permanova_permdisp.tsv`; two rarefied tables: `table-analysis-rarefied.qza` (base) and `{db}-table-rarefied.qza` (per-DB)                                                                                                                 |
| 5 — Diff. Abundance | LEfSe LDA + cladogram plots (`.svg` + `.png`) from the vanilla implementation (`lefse_Huttenhower_*`); ANCOMBC2 effect-size + cladogram pair rendered through the same shared renderer (`ANCOMBC2_EffectSize_by_*`, `ANCOMBC2_Cladogram_by_*`, log2), plus the raw results `.tsv` |
| 6 — Function        | PICRUSt2 KO/EC/pathway TSVs,`picrust2_KO`/`picrust2_EC`/`picrust2_MetaCyc` heatmaps (`.svg` + `.png`)                                                                                                                                                                               |
| 7 — Report          | `visualisations_report.pdf`                                                                                                                                                                                                                                                                 |

> **Note on plot formats:** all plots are emitted as paired `.svg` + `.png` files.
