
rule NGS_import:
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
        "[QIIME] Demultiplexed NGS sequences found. Importing as QIIME Artifact..."
    log:
        STUDY_DIR / "logs/NGS_import.log"
    shell:
        """
        exec > {log} 2>&1
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-format PairedEndFastqManifestPhred33V2 \
            --input-path {input.manifest} \
            --output-path {output.demux_qza}
        """

rule NGS_cutadapt:
    input:
        demux_qza = QIIME_DIR / "demux.qza"
    output:
        trimmed_qza = QIIME_DIR / "trimmed-demux.qza"
    container:
        QIIME_CONTAINER
    threads: n_threads
    resources:
        mem_mb = 8000
    params:
        error_rate = cutadapt_error_rate
    message:
        "[QIIME] Trimming primers using cutadapt..."
    log:
        STUDY_DIR / "logs/NGS_cutadapt.log"
    shell:
        """
        exec > {log} 2>&1
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {input.demux_qza} \
            --p-cores {threads} \
            --p-error-rate {params.error_rate} \
            --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
            --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
            --o-trimmed-sequences {output.trimmed_qza}
        """

rule NGS_summarize_trimmed_demux:
    input:
        trimmed_qza = QIIME_DIR / "trimmed-demux.qza"
    output:
        trimmed_quality_qzv = QIIME_DIR / "trimmed-quality.qzv"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        "[QIIME] Generating summary of trimmed demultiplexed NGS sequences..."
    log:
        STUDY_DIR / "logs/NGS_summarize_trimmed_demux.log"
    shell:
        """
        exec > {log} 2>&1
        qiime demux summarize \
            --i-data {input.trimmed_qza} \
            --o-visualization {output.trimmed_quality_qzv}
        """

rule NGS_export_trimmed_quality:
    input:
        trimmed_quality_qzv = QIIME_DIR / "trimmed-quality.qzv"
    output:
        sentinel = QIIME_DIR / "trimmed-quality-tsv/.export_complete"  # Add this!
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    params:
        quality_tsv_dir = QIIME_DIR / "trimmed-quality-tsv"
    message:
        "[QIIME] Exporting NGS quality summary TSV to determine truncation length..."
    log:
        STUDY_DIR / "logs/NGS_export_trimmed_quality.log"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.quality_tsv_dir}
        qiime tools export \
            --input-path {input.trimmed_quality_qzv} \
            --output-path {params.quality_tsv_dir}
        touch {output.sentinel}
        """

rule NGS_generate_trunc_len:
    input:
        sentinel = QIIME_DIR / "trimmed-quality-tsv/.export_complete"
    output:
        trunc_len_csv = TABLES_DIR / "trunc_len.csv"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    params:
        quality_tsv_dir = QIIME_DIR / "trimmed-quality-tsv",
        q_threshold = trunc_len_q_threshold
    message:
        "[QIIME] Generating NGS truncation length based on quality summary..."
    log:
        STUDY_DIR / "logs/NGS_generate_trunc_len.log"
    shell:
        """
        exec > {log} 2>&1
        python scripts/generate_trunc_length.py \
            {params.quality_tsv_dir}/data.jsonp \
            {output.trunc_len_csv} {params.q_threshold}
        """


def read_trunc_len(wildcards):
    import csv
    csv_file = TABLES_DIR / "trunc_len.csv"
    try:
        with open(csv_file) as f:
            reader = csv.DictReader(f)
            row = next(reader)
            return {"trunc_len_f": int(row["trunc_len_f"]),
                    "trunc_len_r": int(row["trunc_len_r"])}
    except (StopIteration, FileNotFoundError): # used for touching
        return {"trunc_len_f": 0, "trunc_len_r": 0}

# Base names
BASE_TABLE = "table-dada2.qza"
BASE_REP = "rep-seqs-dada2.qza"
DADA2_STATS = "dada2-stats.qza"

BASE_TABLE_QZV = BASE_TABLE.replace(".qza", ".qzv")
BASE_REP_QZV = BASE_REP.replace(".qza", ".qzv")
DADA2_STATS_QZV = DADA2_STATS.replace(".qza", ".qzv")


rule NGS_dada2:
    input:
        trimmed_qza = QIIME_DIR / "trimmed-demux.qza",
        trunc_len_csv = TABLES_DIR / "trunc_len.csv"
    output:
        table = QIIME_DIR / BASE_TABLE,
        rep_seqs = QIIME_DIR / BASE_REP,
        stats = QIIME_DIR / DADA2_STATS,
        base_transition = QIIME_DIR / "base-transition-stats-dada2.qza"
    params:
        trim_left_f = dada2_trim_left_f,
        trim_left_r = dada2_trim_left_r,
        trunc_len_f = lambda wildcards: read_trunc_len(wildcards)["trunc_len_f"],
        trunc_len_r = lambda wildcards: read_trunc_len(wildcards)["trunc_len_r"]
    threads: n_threads
    container:
        QIIME_CONTAINER
    # runtime: mirrors TGS_dada2's reasoning -- scale with retry attempt
    # rather than a fixed guess, since actual denoising time varies by
    # project dataset size (see TGS_dada2 for the observed range).
    resources:
        mem_mb = 16000,
        runtime = lambda wildcards, attempt: 240 * attempt
    message:
        "[QIIME] Running DADA2 denoising for NGS data..."
    log:
        STUDY_DIR / "logs/NGS_dada2.log"
    shell:
        """
        exec > {log} 2>&1
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input.trimmed_qza} \
            --p-trim-left-f {params.trim_left_f} \
            --p-trim-left-r {params.trim_left_r} \
            --p-trunc-len-f {params.trunc_len_f} \
            --p-trunc-len-r {params.trunc_len_r} \
            --o-table {output.table} \
            --o-representative-sequences {output.rep_seqs} \
            --o-denoising-stats {output.stats} \
            --o-base-transition-stats {output.base_transition} \
            --p-n-threads {threads}
        """



rule NGS_visualise_dada2_outputs:
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
        STUDY_DIR / "logs/NGS_visualise_dada2_outputs.log"
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

rule NGS_export_table_summary:
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
        "[QIIME] Exporting NGS table summary..."
    log:
        STUDY_DIR / "logs/NGS_export_table_summary.log"
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


rule NGS_generate_rarefy_depth:
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
        "[PYTHON] Reading NGS non-phylogenetic rarefaction depth from QIIME summary..."
    log:
        STUDY_DIR / "logs/NGS_generate_rarefy_depth.log"
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

rule NGS_make_rarefied_version:
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
        "[QIIME] Generating NGS rarefied feature table for non-phylogenetic diversity metrics..."
    log:
        STUDY_DIR / "logs/NGS_make_rarefied_version.log"
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


rule NGS_export_table_to_biom:
    input:
        table_qza = QIIME_DIR / "table-analysis.qza"
    output:
        table_biom = TABLES_DIR / "study-seqs.biom"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[QIIME] Exporting NGS feature table to BIOM format..."
    log:
        STUDY_DIR / "logs/NGS_export_table_to_biom.log"
    shell:
        """
        exec > {log} 2>&1
        qiime tools export \
            --input-path {input.table_qza} \
            --output-path exported_table_temp
        mv exported_table_temp/feature-table.biom {output.table_biom}
        rm -r exported_table_temp
        """

rule NGS_export_rep_seqs_to_fna:
    input:
        rep_seqs_qza = QIIME_DIR / "rep-seqs-analysis.qza"
    output:
        rep_seqs_fna = TABLES_DIR / "study-seqs.fna"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[QIIME] Exporting NGS representative sequences to FASTA..."
    log:
        STUDY_DIR / "logs/NGS_export_rep_seqs_to_fna.log"
    shell:
        """
        exec > {log} 2>&1
        qiime tools export \
            --input-path {input.rep_seqs_qza} \
            --output-path exported_seqs_temp
        mv exported_seqs_temp/dna-sequences.fasta {output.rep_seqs_fna}
        rm -r exported_seqs_temp
        """
