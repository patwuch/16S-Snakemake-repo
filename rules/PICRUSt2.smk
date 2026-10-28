
PICRUST2_CONTAINER = str(WORKFLOW_DIR / "containers/picrust2-env.sif")

rule picrust2_run:
    input:
        seqs_fna = TABLES_DIR / "study-seqs.fna",
        table_biom = TABLES_DIR / "study-seqs.biom"
    output:
        picrust2_dir = directory(STUDY_DIR / "picrust2_output")
    threads: 12
    container: PICRUST2_CONTAINER
    resources:
        mem_mb = 16000
    message:
        """
        [PICRUSt2] Running PICRUSt2 to predict functional profiles (ECs, KOs, Pathways) from ASV-level FASTA and BIOM...
        """
    log:
        STUDY_DIR / "logs/picrust2_run.log"
    shell:
        """
        exec > {log} 2>&1
        picrust2_pipeline.py \
            -s {input.seqs_fna} \
            -i {input.table_biom} \
            -o {output.picrust2_dir} \
            -p {threads}
        """

rule picrust2_describe:
    input:
        picrust2_dir = STUDY_DIR / "picrust2_output"
    output:
        ko_described = STUDY_DIR / "picrust2_described" / "KO_metagenome_unstrat_described.tsv.gz",
        ec_described = STUDY_DIR / "picrust2_described" / "EC_metagenome_unstrat_described.tsv.gz",
        pathway_described = STUDY_DIR / "picrust2_described" / "pathway_abun_unstrat_described.tsv.gz"
    threads: 12
    container: PICRUST2_CONTAINER
    resources:
        mem_mb = 4000
    message:
        """
        [PICRUSt2] Adding descriptions to PICRUSt2 predicted functional profiles (ECs, KOs, Pathways)...
        """
    log:
        STUDY_DIR / "logs/picrust2_describe.log"
    shell:
        """
        exec > {log} 2>&1
        mkdir -p $(dirname {output.ko_described})

        add_descriptions.py \
            -i {input.picrust2_dir}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
            -o {output.ko_described} \
            -m KO

        add_descriptions.py \
            -i {input.picrust2_dir}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
            -o {output.ec_described} \
            -m EC

        add_descriptions.py \
            -i {input.picrust2_dir}/pathways_out/path_abun_unstrat.tsv.gz \
            -o {output.pathway_described} \
            -m METACYC
        """

rule picrust2_plot:
    input:
        ec = f"{STUDY_DIR}/picrust2_described/EC_metagenome_unstrat_described.tsv.gz",
        ko = f"{STUDY_DIR}/picrust2_described/KO_metagenome_unstrat_described.tsv.gz",
        pathway = f"{STUDY_DIR}/picrust2_described/pathway_abun_unstrat_described.tsv.gz",
        metadata = f"{STUDY_DIR}/metadata.tsv"
    output:
        ko_svg = f"{STUDY_DIR}/plots/picrust2_KO.svg",
        ec_svg = f"{STUDY_DIR}/plots/picrust2_EC.svg",
        metacyc_svg = f"{STUDY_DIR}/plots/picrust2_MetaCyc.svg",
        ko_png = f"{STUDY_DIR}/plots/picrust2_KO.png",
        ec_png = f"{STUDY_DIR}/plots/picrust2_EC.png",
        metacyc_png = f"{STUDY_DIR}/plots/picrust2_MetaCyc.png"
    params:
        top_n = top_n_picrust
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        """
        [PICRUSt2] Plotting PICRUSt2 predicted functional profiles (ECs, KOs, Pathways) as a heatmap for the top {params.top_n} most abundant features across samples..."""
    log:
        STUDY_DIR / "logs/picrust2_plot.log"
    script:
        "../scripts/plot_picrust2.R"