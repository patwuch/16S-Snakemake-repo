rule TGS_import:
    input:
        manifest = STUDY_DIR / "manifest.tsv",
        make_dirs_marker = STUDY_DIR / ".dirs_created"
    output:
        demux_qza = QIIME_DIR / "demux.qza"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[QIIME] Demultiplexed TGS sequences found. Importing as QIIME Artifact..."
    log:
        STUDY_DIR / "logs/TGS_import.log"
    shell:
        """
        exec > {log} 2>&1
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-format SingleEndFastqManifestPhred33V2 \
            --input-path {input.manifest} \
            --output-path {output.demux_qza}
        echo "[QIIME] TGS Circular consensus sequences leaves no adapter. Skipping cutadapt..."
        echo "[INFO] TGS Circular consensus sequences do not need to be truncated, setting trim length to 0."
        """

rule TGS_summarize_demux:
    input:
        demux_qza = QIIME_DIR / "demux.qza"
    output:
        demux_quality_qzv = QIIME_DIR / "demux-quality.qzv"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        "[QIIME] Generating summary of demultiplexed TGS sequences..."
    log:
        STUDY_DIR / "logs/TGS_summarize_demux.log"
    shell:
        """
        exec > {log} 2>&1
        qiime demux summarize \
            --i-data {input.demux_qza} \
            --o-visualization {output.demux_quality_qzv}
        """

# Base names
BASE_TABLE = "table-dada2.qza"
BASE_REP = "rep-seqs-dada2.qza"
DADA2_STATS = "dada2-stats.qza"

BASE_TABLE_QZV = BASE_TABLE.replace(".qza", ".qzv")
BASE_REP_QZV = BASE_REP.replace(".qza", ".qzv")
DADA2_STATS_QZV = DADA2_STATS.replace(".qza", ".qzv")


# Circular consensus does not need to be truncated so set length to 0
rule TGS_dada2:
    input:
        demux_qza = QIIME_DIR / "demux.qza"
    output:
        table = QIIME_DIR / BASE_TABLE,
        repseqs = QIIME_DIR / BASE_REP,
        stats = QIIME_DIR / DADA2_STATS,
        base_transition = QIIME_DIR / "base-transition-stats-dada2.qza"
    threads: n_threads
    container:
        QIIME_CONTAINER
    # runtime: observed 44min (tisha260707) to >60min (hu260708, timed out
    # under the old shared 60min profile default) for CCS denoising -- real
    # per-project variance, so scale with retry attempt rather than picking
    # one fixed number that has to cover every future project's dataset size.
    resources:
        mem_mb = 16000,
        runtime = lambda wildcards, attempt: 240 * attempt
    message:
        "[QIIME] Running DADA2 denoising for TGS data..."
    log:
        STUDY_DIR / "logs/TGS_dada2.log"
    shell:
        """
        exec > {log} 2>&1
        qiime dada2 denoise-ccs \
            --i-demultiplexed-seqs {input.demux_qza} \
            --p-front AGAGTTTGATCMTGGCTCAG \
            --o-table {output.table} \
            --o-representative-sequences {output.repseqs} \
            --o-denoising-stats {output.stats} \
            --o-base-transition-stats {output.base_transition} \
            --p-n-threads {threads} \
            --p-trunc-len 0
        """


rule TGS_visualise_dada2_outputs:
    input:
        table = QIIME_DIR / BASE_TABLE,
        repseqs = QIIME_DIR / BASE_REP,
        stats = QIIME_DIR / DADA2_STATS
    output:
        table_qzv = QIIME_DIR / BASE_TABLE_QZV,
        repseqs_qzv = QIIME_DIR / BASE_REP_QZV,
        stats_qzv = QIIME_DIR / DADA2_STATS_QZV
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        "[QIIME] Generating visualizations of DADA2 outputs..."
    log:
        STUDY_DIR / "logs/TGS_visualise_dada2_outputs.log"
    shell:
        """
        exec > {log} 2>&1
        qiime feature-table summarize \
            --i-table {input.table} \
            --o-visualization {output.table_qzv}

        qiime feature-table tabulate-seqs \
            --i-data {input.repseqs} \
            --o-visualization {output.repseqs_qzv}

        qiime metadata tabulate \
            --m-input-file {input.stats} \
            --o-visualization {output.stats_qzv}
        """

rule TGS_export_table_summary:
    input:
        table_qza = QIIME_DIR / "table-analysis.qza"
    output:
        summary_tsv = QIIME_DIR / "table-summary/feature-table.tsv",
        sentinel = QIIME_DIR / "table-summary/.export_complete"
    params:
        outdir = QIIME_DIR / "table-summary"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[QIIME] Exporting TGS table summary..."
    log:
        STUDY_DIR / "logs/TGS_export_table_summary.log"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        qiime tools export \
            --input-path {input.table_qza} \
            --output-path {params.outdir}
        biom convert \
            -i {params.outdir}/feature-table.biom \
            -o {output.summary_tsv} \
            --to-tsv
        touch {output.sentinel}
        """

rule TGS_generate_rarefy_depth_nonphylogenetic:
    input:
        sentinel = QIIME_DIR / "table-summary/.export_complete",
        summary_tsv = QIIME_DIR / "table-summary/feature-table.tsv"
    output:
        rarefy_csv = TABLES_DIR / "rarefy_depth.csv"
    params:
        percentile = rarefy_depth_percentile
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[PYTHON] Reading TGS non-phylogenetic rarefaction depth from QIIME summary..."
    log:
        STUDY_DIR / "logs/TGS_generate_rarefy_depth_nonphylogenetic.log"
    shell:
        """
        exec > {log} 2>&1
        python - << EOF
import pandas as pd
import numpy as np
summary_file = "{input.summary_tsv}"
out_file = "{output.rarefy_csv}"
percentile = {params.percentile}

df = pd.read_csv(summary_file, skiprows=1, sep='\t', index_col=0)
sample_sums = df.sum(axis=0)  # sum across features for each sample
if percentile == 0:
    depth = int(sample_sums.min())
else:
    depth = int(np.percentile(sample_sums, percentile))
pd.DataFrame([{{"rarefaction_depth": depth}}]).to_csv(out_file, index=False)
EOF
        """

rule TGS_generate_rarefied_version_nonphylogenetic:
    input:
        table = QIIME_DIR / "table-analysis.qza",
        rarefy_csv = TABLES_DIR / "rarefy_depth.csv"
    output:
        rarefied_table = QIIME_DIR / "table-analysis-rarefied.qza"
    params:
        depth = lambda wildcards: read_rarefy_depth(wildcards)
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        "[QIIME] Generating TGS rarefied feature table for non-phylogenetic diversity metrics..."
    log:
        STUDY_DIR / "logs/TGS_generate_rarefied_version_nonphylogenetic.log"
    shell:
        """
        exec > {log} 2>&1
        qiime feature-table rarefy \
            --i-table {input.table} \
            --p-sampling-depth {params.depth} \
            --o-rarefied-table {output.rarefied_table}
        """

def read_rarefy_depth(wildcards):
    import csv
    csv_file = TABLES_DIR / "rarefy_depth.csv"
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        row = next(reader)
        return int(row["rarefaction_depth"])


rule TGS_export_table_to_biom:
    input:
        table_qza = QIIME_DIR / "table-analysis.qza"
    output:
        table_biom = TABLES_DIR / "study-seqs.biom"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[QIIME] Exporting TGS feature table to BIOM format..."
    log:
        STUDY_DIR / "logs/TGS_export_table_to_biom.log"
    shell:
        """
        exec > {log} 2>&1
        qiime tools export \
            --input-path {input.table_qza} \
            --output-path exported_table_temp
        mv exported_table_temp/feature-table.biom {output.table_biom}
        rm -r exported_table_temp
        """

rule TGS_export_rep_seqs_to_fna:
    input:
        rep_seqs_qza = QIIME_DIR / "rep-seqs-analysis.qza"
    output:
        rep_seqs_fna = TABLES_DIR / "study-seqs.fna"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[QIIME] Exporting TGS representative sequences to FASTA..."
    log:
        STUDY_DIR / "logs/TGS_export_rep_seqs_to_fna.log"
    shell:
        """
        exec > {log} 2>&1
        qiime tools export \
            --input-path {input.rep_seqs_qza} \
            --output-path exported_seqs_temp
        mv exported_seqs_temp/dna-sequences.fasta {output.rep_seqs_fna}
        rm -r exported_seqs_temp
        """
