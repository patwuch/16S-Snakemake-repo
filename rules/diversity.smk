
rule calc_alpha_diversity_non_phylogenetic:
    input:
        table = QIIME_DIR / "table-analysis-rarefied.qza"
    output:
        shannon  = str(CORE_METRICS_DIR / "shannon-vector.qza"),
        simpson  = str(CORE_METRICS_DIR / "simpson-vector.qza"),
        evenness = str(CORE_METRICS_DIR / "evenness-vector.qza"),
        chao1    = str(CORE_METRICS_DIR / "chao1-vector.qza")
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        """
        [QIIME] Calculating non-phylogenetic alpha diversity metrics (Shannon, Simpson, Pielou's Evenness, Chao1) from rarefied feature table...
        """
    log:
        STUDY_DIR / "logs/calc_alpha_diversity_non_phylogenetic.log"
    shell:
        """
        exec > {log} 2>&1
        qiime diversity alpha --i-table {input.table} --p-metric shannon --o-alpha-diversity {output.shannon}
        qiime diversity alpha --i-table {input.table} --p-metric simpson --o-alpha-diversity {output.simpson}
        qiime diversity alpha --i-table {input.table} --p-metric pielou_e --o-alpha-diversity {output.evenness}
        qiime diversity alpha --i-table {input.table} --p-metric chao1 --o-alpha-diversity {output.chao1}
        """

rule calc_beta_diversity_non_phylogenetic:
    input:
        table_rarefied = QIIME_DIR / "table-analysis-rarefied.qza",
    output:
        jaccard    = str(CORE_METRICS_DIR / "jaccard-distance-matrix.qza"),
        braycurtis = str(CORE_METRICS_DIR / "bray-curtis-distance-matrix.qza"),
    container:
        QIIME_CONTAINER
    threads: n_threads
    resources:
        mem_mb = 4000
    message:
        """
        [QIIME] Calculating non-phylogenetic beta diversity distance matrices (Jaccard, Bray-Curtis) from rarefied feature table..."""
    log:
        STUDY_DIR / "logs/calc_beta_diversity_non_phylogenetic.log"
    shell:
        """
        exec > {log} 2>&1
        qiime diversity beta \
            --i-table {input.table_rarefied} \
            --p-metric jaccard \
            --p-n-jobs {threads} \
            --o-distance-matrix {output.jaccard}
        qiime diversity beta \
            --i-table {input.table_rarefied} \
            --p-metric braycurtis \
            --p-n-jobs {threads} \
            --o-distance-matrix {output.braycurtis}
        """



def find_phylogeny(wildcards):
    # Silva138 has no static backbone tree shipped with the reference database —
    # its tree is built per-study by fragment-insertion (SEPP) in Silva138.smk,
    # so it lives under QIIME_DIR rather than the shared REF_DIR.
    if wildcards.db == "Silva138":
        return str(QIIME_DIR / "Silva138-tree.qza")
    pattern = os.path.join(REF_DIR, wildcards.db, "*.nwk.qza")
    matches = glob.glob(pattern)
    if len(matches) == 0:
        raise ValueError(f"No .nwk.qza phylogeny found for db: {wildcards.db}")
    if len(matches) > 1:
        raise ValueError(f"Multiple .nwk.qza phylogenies found for db: {wildcards.db}: {matches}")
    return matches[0]



def read_db_rarefy_depth(wildcards):
    import csv
    csv_file = TABLES_DIR / f"{wildcards.db}_rarefy_depth.csv"
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        row = next(reader)
        return int(row["rarefaction_depth"])


rule export_db_table_for_summary_tsv:
    input:
        table = QIIME_DIR / "{db}-table.qza"
    output:
        summary_tsv = temp(TABLES_DIR / "{db}-table-summary" / "feature-table.tsv")
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        """
        [QIIME] Exporting {wildcards.db} feature table to TSV for summary statistics calculation...
        """
    log:
        STUDY_DIR / "logs/export_db_table_for_summary_tsv/{db}.log"
    shell:
        """
        exec > {log} 2>&1
        tmpdir=$(mktemp -d)
        qiime tools export \
            --input-path {input.table} \
            --output-path $tmpdir
        biom convert \
            -i $tmpdir/feature-table.biom \
            -o {output.summary_tsv} \
            --to-tsv
        rm -rf $tmpdir
        """


rule generate_db_rarefy_depth:
    input:
        summary_tsv = TABLES_DIR / "{db}-table-summary" / "feature-table.tsv"
    output:
        rarefy_csv = TABLES_DIR / "{db}_rarefy_depth.csv"
    params:
        percentile = rarefy_depth_percentile
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        """
        [QIIME] Calculating rarefaction depth for {wildcards.db} feature table based on {params.percentile} percentile of sample sums...
        """
    log:
        STUDY_DIR / "logs/generate_db_rarefy_depth/{db}.log"
    shell:
        """
        exec > {log} 2>&1
        python - << 'EOF'
import pandas as pd
import numpy as np

df = pd.read_csv("{input.summary_tsv}", skiprows=1, sep='\t', index_col=0)
sample_sums = df.sum(axis=0)
percentile = {params.percentile}
depth = int(sample_sums.min()) if percentile == 0 else int(np.percentile(sample_sums, percentile))
pd.DataFrame([{{"rarefaction_depth": depth}}]).to_csv("{output.rarefy_csv}", index=False)
EOF
        """


rule rarefy_phylogenetic_table:
    input:
        table      = QIIME_DIR / "{db}-table.qza",
        rarefy_csv = TABLES_DIR / "{db}_rarefy_depth.csv"
    output:
        rarefied = QIIME_DIR / "{db}-table-rarefied.qza"
    params:
        depth = lambda wildcards: read_db_rarefy_depth(wildcards)
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        """
        [QIIME] Rarefying {wildcards.db} feature table to depth {params.depth}...
        """
    log:
        STUDY_DIR / "logs/rarefy_phylogenetic_table/{db}.log"
    shell:
        """
        exec > {log} 2>&1
        qiime feature-table rarefy \
            --i-table {input.table} \
            --p-sampling-depth {params.depth} \
            --o-rarefied-table {output.rarefied}
        """

rule calc_alpha_diversity_phylogenetic:
    input:
        table     = QIIME_DIR / "{db}-table-rarefied.qza",
        phylogeny = find_phylogeny
    output:
        faith_pd = str(CORE_METRICS_DIR / "{db}-faith-pd-vector.qza")
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        """
        [QIIME] Calculating phylogenetic alpha diversity metric (Faith's PD) from rarefied phylogenetic feature table...
        """
    log:
        STUDY_DIR / "logs/calc_alpha_diversity_phylogenetic/{db}.log"
    shell:
        """
        exec > {log} 2>&1
        qiime diversity alpha-phylogenetic \
            --i-table {input.table} \
            --i-phylogeny {input.phylogeny} \
            --p-metric faith_pd \
            --o-alpha-diversity {output.faith_pd}
        """


rule calc_beta_diversity_phylogenetic:
    input:
        table_rarefied = QIIME_DIR / "{db}-table-rarefied.qza",
        phylogeny      = find_phylogeny
    output:
        unweighted = str(CORE_METRICS_DIR / "{db}-unweighted-unifrac-distance-matrix.qza"),
        weighted   = str(CORE_METRICS_DIR / "{db}-weighted-unifrac-distance-matrix.qza")
    container:
        QIIME_CONTAINER
    threads: n_threads
    resources:
        mem_mb = 8000
    message:
        """
        [QIIME] Calculating phylogenetic beta diversity distance matrices (Unweighted and Weighted UniFrac) from rarefied phylogenetic feature table...
        """
    log:
        STUDY_DIR / "logs/calc_beta_diversity_phylogenetic/{db}.log"
    shell:
        """
        exec > {log} 2>&1
        qiime diversity beta-phylogenetic \
            --i-table {input.table_rarefied} \
            --i-phylogeny {input.phylogeny} \
            --p-metric unweighted_unifrac \
            --p-threads {threads} \
            --o-distance-matrix {output.unweighted}
        qiime diversity beta-phylogenetic \
            --i-table {input.table_rarefied} \
            --i-phylogeny {input.phylogeny} \
            --p-metric weighted_unifrac \
            --p-threads {threads} \
            --o-distance-matrix {output.weighted}
        """
rule pcoa_beta_diversity_non_phylogenetic:
    input:
        jaccard    = str(CORE_METRICS_DIR / "jaccard-distance-matrix.qza"),
        braycurtis = str(CORE_METRICS_DIR / "bray-curtis-distance-matrix.qza"),
    output:
        jaccard_pcoa    = str(CORE_METRICS_DIR / "jaccard-pcoa-results.qza"),
        braycurtis_pcoa = str(CORE_METRICS_DIR / "bray-curtis-pcoa-results.qza"),
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        """
        [QIIME] Calculating non-phylogenetic beta diversity distance matrices (Jaccard, Bray-Curtis) from distance matrices..."""
    log:
        STUDY_DIR / "logs/pcoa_beta_diversity_non_phylogenetic.log"
    shell:
        """
        exec > {log} 2>&1
        qiime diversity pcoa \
            --i-distance-matrix {input.braycurtis} \
            --o-pcoa {output.braycurtis_pcoa}
        qiime diversity pcoa \
            --i-distance-matrix {input.jaccard} \
            --o-pcoa {output.jaccard_pcoa}
        """

rule pcoa_beta_diversity_phylogenetic:
    input:
        unweighted = str(CORE_METRICS_DIR / "{db}-unweighted-unifrac-distance-matrix.qza"),
        weighted   = str(CORE_METRICS_DIR / "{db}-weighted-unifrac-distance-matrix.qza"),
    output:
        unweighted_pcoa = str(CORE_METRICS_DIR / "{db}-unweighted-unifrac-pcoa-results.qza"),
        weighted_pcoa   = str(CORE_METRICS_DIR / "{db}-weighted-unifrac-pcoa-results.qza"),
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        """
        [QIIME] Calculating phylogenetic beta diversity distance matrices (Unweighted and Weighted UniFrac) from distance matrices...
        """
    log:
        STUDY_DIR / "logs/pcoa_beta_diversity_phylogenetic/{db}.log"
    shell:
        """
        exec > {log} 2>&1
        qiime diversity pcoa \
            --i-distance-matrix {input.unweighted} \
            --o-pcoa {output.unweighted_pcoa}
        qiime diversity pcoa \
            --i-distance-matrix {input.weighted} \
            --o-pcoa {output.weighted_pcoa}
        """

# Standalone per-metric alpha diversity figures. Shannon/Evenness/Simpson/Chao1
# don't need a tree; Faith PD does, so it's only ever requested (via
# alpha_phylogenetic_outputs / alpha_metric_outputs in the Snakefile) for dbs
# in phylogenetic_reference_db.
ALPHA_METRIC_VECTOR = {
    "shannon":  lambda db: CORE_METRICS_DIR / "shannon-vector.qza",
    "evenness": lambda db: CORE_METRICS_DIR / "evenness-vector.qza",
    "simpson":  lambda db: CORE_METRICS_DIR / "simpson-vector.qza",
    "chao1":    lambda db: CORE_METRICS_DIR / "chao1-vector.qza",
    "faith_pd": lambda db: CORE_METRICS_DIR / f"{db}-faith-pd-vector.qza",
}

def alpha_metric_vector_input(wildcards):
    return str(ALPHA_METRIC_VECTOR[wildcards.metric](wildcards.db))

rule plot_alpha_diversity_metric:
    wildcard_constraints:
        metric = "shannon|evenness|simpson|chao1|faith_pd"
    input:
        vector   = alpha_metric_vector_input,
        metadata = STUDY_DIR / "metadata.tsv"
    output:
        plot_svg = ALPHA_DIR / "{db}_{metric}_{group_col}.svg",
        plot_png = ALPHA_DIR / "{db}_{metric}_{group_col}.png",
        sentinel = DIVERSITY_DIR / ".{db}_{metric}_{group_col}_done"
    params:
        group_by      = "{group_col}",
        db            = "{db}",
        metric        = "{metric}",
        group_order   = lambda wc: GROUP_ORDERS[wc.group_col],
        color_palette = lambda wc: GROUP_COLORS[wc.group_col]
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        """
        [QIIME] Plotting {wildcards.metric} alpha diversity for {wildcards.db} grouped by {params.group_by}...
        """
    log:
        STUDY_DIR / "logs/plot_alpha_diversity_metric/{db}_{group_col}_{metric}.log"
    script:
        SCRIPTS_DIR / "plot_alpha_diversity_metric.py"

# Standalone per-metric beta diversity figures. Jaccard/Bray-Curtis don't need
# a tree; the two UniFrac metrics do, so they're only ever requested for dbs
# in phylogenetic_reference_db.
BETA_METRIC_PCOA = {
    "jaccard":            lambda db: CORE_METRICS_DIR / "jaccard-pcoa-results.qza",
    "braycurtis":         lambda db: CORE_METRICS_DIR / "bray-curtis-pcoa-results.qza",
    "weighted_unifrac":   lambda db: CORE_METRICS_DIR / f"{db}-weighted-unifrac-pcoa-results.qza",
    "unweighted_unifrac": lambda db: CORE_METRICS_DIR / f"{db}-unweighted-unifrac-pcoa-results.qza",
}

def beta_metric_pcoa_input(wildcards):
    return str(BETA_METRIC_PCOA[wildcards.metric](wildcards.db))

rule plot_beta_diversity_metric:
    wildcard_constraints:
        metric = "jaccard|braycurtis|weighted_unifrac|unweighted_unifrac"
    input:
        pcoa     = beta_metric_pcoa_input,
        metadata = STUDY_DIR / "metadata.tsv"
    output:
        plot_svg = BETA_DIR / "{db}_beta_{metric}_{group_col}.svg",
        plot_png = BETA_DIR / "{db}_beta_{metric}_{group_col}.png",
        sentinel = DIVERSITY_DIR / ".{db}_beta_{metric}_{group_col}_done"
    params:
        group_by      = "{group_col}",
        db            = "{db}",
        metric        = "{metric}",
        group_order   = lambda wc: GROUP_ORDERS[wc.group_col],
        color_palette = lambda wc: GROUP_COLORS[wc.group_col]
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        """
        [QIIME] Plotting {wildcards.metric} beta diversity PCoA for {wildcards.db} grouped by {params.group_by}...
        """
    log:
        STUDY_DIR / "logs/plot_beta_diversity_metric/{db}_{group_col}_{metric}.log"
    script:
        SCRIPTS_DIR / "plot_beta_diversity_metric.py"

def permanova_dist_inputs(wildcards):
    dist = [
        str(CORE_METRICS_DIR / "jaccard-distance-matrix.qza"),
        str(CORE_METRICS_DIR / "bray-curtis-distance-matrix.qza"),
    ]
    if wildcards.db in phylogenetic_reference_db:
        dist += [
            str(CORE_METRICS_DIR / f"{wildcards.db}-unweighted-unifrac-distance-matrix.qza"),
            str(CORE_METRICS_DIR / f"{wildcards.db}-weighted-unifrac-distance-matrix.qza"),
        ]
    return dist

rule run_permanova_permdisp:
    input:
        dist = permanova_dist_inputs,
        meta = STUDY_DIR / "metadata.tsv"
    output:
        TABLES_DIR / "{db}_permanova_permdisp.tsv"
    container:
        QIIME_CONTAINER
    threads: n_threads
    resources:
        mem_mb = 4000
    message:
        """
        [QIIME] Running PERMANOVA and PERMDISP tests for {wildcards.db} beta diversity distance matrices...
        """
    log:
        STUDY_DIR / "logs/run_permanova_permdisp/{db}.log"
    script:
        str(SCRIPTS_DIR / "run_permanova_permdisp.py")