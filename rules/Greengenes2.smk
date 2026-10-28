GREENGENES2_CONTAINER = str(WORKFLOW_DIR / "containers/qiime2-2025.10-amplicon-Greengenes2.sif")


if region == "region_V3V4":
    rule GG2_V3V4_taxonomy:
        input:
            rep_seqs = QIIME_DIR / "rep-seqs-analysis.qza",
            greengenes_nb_classifier = REF_DIR / "Greengenes2" / "2024.09.custom.V3V4.nb.qza"
        output:
            taxonomy = QIIME_DIR / "Greengenes2-taxonomy.qza",
            taxonomy_qzv = QIIME_DIR / "Greengenes2-taxonomy.qzv",
            sentinel = QIIME_DIR / ".Greengenes2_taxonomy_done"
        container:
            GREENGENES2_CONTAINER
        resources:
            mem_mb = 20000
        message:
            "[GG2] Classifying 16S region V3V4 sequences using Greengenes2 Naive Bayes Classifier..."
        log:
            STUDY_DIR / "logs/GG2_V3V4_taxonomy.log"
        shell:
            """
            exec > {log} 2>&1
            qiime feature-classifier classify-sklearn \
                --i-classifier {input.greengenes_nb_classifier} \
                --i-reads {input.rep_seqs} \
                --o-classification {output.taxonomy}
            echo "[GG2] Generating taxonomy visualization..."
            qiime metadata tabulate \
                --m-input-file {output.taxonomy} \
                --o-visualization {output.taxonomy_qzv}
            touch {output.sentinel}
            """

elif region == "full_length":
    rule GG2_Full_taxonomy:
        input:
            rep_seqs = QIIME_DIR / "rep-seqs-analysis.qza",
            greengenes_nb_classifier = REF_DIR / "Greengenes2" / "2024.09.backbone.full-length.nb.qza"
        output:
            taxonomy = QIIME_DIR / "Greengenes2-taxonomy.qza",
            taxonomy_qzv = QIIME_DIR / "Greengenes2-taxonomy.qzv",
            sentinel = QIIME_DIR / ".Greengenes2_taxonomy_done"
        container:
            GREENGENES2_CONTAINER
        resources:
            mem_mb = 20000
        message:
            "[GG2] Classifying 16S Full-length sequences using Greengenes2 Naive Bayes Classifier..."
        log:
            STUDY_DIR / "logs/GG2_Full_taxonomy.log"
        shell:
            """
            exec > {log} 2>&1
            qiime feature-classifier classify-sklearn \
                --i-classifier {input.greengenes_nb_classifier} \
                --i-reads {input.rep_seqs} \
                --o-classification {output.taxonomy}
            echo "[GG2] Generating taxonomy visualization..."
            qiime metadata tabulate \
                --m-input-file {output.taxonomy} \
                --o-visualization {output.taxonomy_qzv}
            touch {output.sentinel}
            """


rule GG2_phylogeny:
    input:
        table = QIIME_DIR / "table-analysis.qza",
        rep_seqs = QIIME_DIR / "rep-seqs-analysis.qza",
        greengenes_db = REF_DIR / "Greengenes2" / "2024.09.backbone.full-length.fna.qza",
        sentinel = QIIME_DIR / ".Greengenes2_taxonomy_done"
    output:
        phylogeny_table = QIIME_DIR / "Greengenes2-table.qza",
        phylogeny_rep_seqs = QIIME_DIR / "Greengenes2-rep-seqs.qza",
        sentinel = QIIME_DIR / ".Greengenes2_phylogeny_done"
    container:
        GREENGENES2_CONTAINER
    resources:
        mem_mb = 16000
    message:
        "[GG2] Building phylogenetic tree for 16S (V3V4/Full-length) sequences using default Greengenes2 SEPP method..."
    log:
        STUDY_DIR / "logs/GG2_phylogeny.log"
    shell:
        """
        exec > {log} 2>&1
        qiime greengenes2 non-v4-16s \
            --i-table {input.table} \
            --i-sequences {input.rep_seqs} \
            --i-backbone {input.greengenes_db} \
            --o-mapped-table {output.phylogeny_table} \
            --o-representatives {output.phylogeny_rep_seqs}
        touch {output.sentinel}
        """

rule GG2_generate_taxa_barplot_qiime:
    input:
        table = QIIME_DIR / "table-analysis.qza",
        taxonomy = QIIME_DIR / "Greengenes2-taxonomy.qza",
        metadata = metadata_path
    output:
        barplot_qzv = QIIME_DIR / "Greengenes2-taxa-bar-plots.qzv",
        sentinel = QIIME_DIR / ".Greengenes2_taxa_barplot_done"
    container:
        GREENGENES2_CONTAINER
    resources:
        mem_mb = 4000
    message:
        "[GG2] Generating Greengenes2 QIIME taxa barplot visualization..."
    log:
        STUDY_DIR / "logs/GG2_generate_taxa_barplot_qiime.log"
    shell:
        """
        exec > {log} 2>&1
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.barplot_qzv}
        touch {output.sentinel}
        """


rule GG2_export_taxonomy_as_tsv:
    input:
        taxonomy_qza = QIIME_DIR / "Greengenes2-taxonomy.qza"
    output:
        taxonomy = TABLES_DIR / "exported-taxonomy" / "Greengenes2_taxonomy.tsv"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[GG2] Exporting Greengenes2 taxonomy as TSV..."
    log:
        STUDY_DIR / "logs/GG2_export_taxonomy_as_tsv.log"
    shell:
        """
        exec > {log} 2>&1
        qiime tools export \
            --input-path {input.taxonomy_qza} \
            --output-path {TABLES_DIR}/exported-taxonomy-temp-gg2

        # Rename exported taxonomy file to your expected filename
        mv {TABLES_DIR}/exported-taxonomy-temp-gg2/taxonomy.tsv \
           {output.taxonomy}
        rm -r {TABLES_DIR}/exported-taxonomy-temp-gg2
        """

# NOTE: the generic taxa-barplot plotting rule (wildcarded on {db}) now lives in
# rules/taxa_barplots.smk, which is included unconditionally — it was moved out
# of this file because Greengenes2.smk is only included when 'Greengenes2' is in
# REFERENCE_DB, which would leave the rule missing for other DBs (e.g. Silva138).