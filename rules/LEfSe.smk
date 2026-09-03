LEFSE_MICROBIOMEMARKER_R_CONTAINER = str(WORKFLOW_DIR / "containers/lefse_MicrobiomeMarker-r-env.sif")

# Both methods end in a LEfSe .res file and share one renderer, so every _plot
# rule runs the same Python script in the same container. Only the _run rules
# still need their method's own R environment.
LEFSE_PLOT_CONTAINER = LEFSE_MICROBIOMEMARKER_R_CONTAINER

# ------------------------------
# Rules
# ------------------------------
rule lefse_MicrobiomeMarker_run:
    input:
        feature_table = DA_FEATURE_TABLE,
        metadata      = STUDY_DIR / "metadata.tsv",
        taxonomy      = TABLES_DIR / "exported-taxonomy" / DA_TAXONOMY_SUFFIX
    output:
        lefse_MicrobiomeMarker_rds = TMP_DIR / "{db}/lefse_MicrobiomeMarker_results_by_{group_col}.rds",
        lefse_MicrobiomeMarker_res = TMP_DIR / "{db}/lefse_MicrobiomeMarker_results_by_{group_col}.res"
    params:
        group_col        = "{group_col}",
        lda_cutoff       = lefse_lda_cutoff,
        kw_cutoff        = lefse_kw_cutoff,
        wilcoxon_cutoff  = lefse_wilcoxon_cutoff,
        norm             = lefse_norm,
        random_seed      = lefse_random_seed
    container:
        LEFSE_MICROBIOMEMARKER_R_CONTAINER
    resources:
        mem_mb = 4000
    script:
        f"{SCRIPTS_DIR}/run_lefse_MicrobiomeMarker.R"

rule lefse_MicrobiomeMarker_plot:
    input:
        lefse_res = TMP_DIR / "{db}/lefse_MicrobiomeMarker_results_by_{group_col}.res"
    output:
        lda_svg       = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_MicrobiomeMarker_LDA_by_{group_col}.svg",
        cladogram_svg = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_MicrobiomeMarker_Cladogram_by_{group_col}.svg",
        lda_png       = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_MicrobiomeMarker_LDA_by_{group_col}.png",
        cladogram_png = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_MicrobiomeMarker_Cladogram_by_{group_col}.png"
    params:
        colors          = lambda wc: GROUP_COLORS[wc.group_col],
        lda_title       = "LEfSe LDA Effect Size",
        cladogram_title = "LEfSe Cladogram"
    container:
        LEFSE_PLOT_CONTAINER
    resources:
        mem_mb = 4000
    script:
        f"{SCRIPTS_DIR}/plot_lefse_res.py"

rule lefse_Huttenhower_run:
    input:
        feature_table = DA_QIIME_TABLE,
        metadata      = STUDY_DIR / "metadata.tsv",
        taxonomy      = DA_QIIME_TAXA
    output:
        lefse_Huttenhower_results = TMP_DIR / "{db}/lefse_Huttenhower_results_by_{group_col}.txt"
    params:
        group_col       = "{group_col}",
        lda_cutoff      = lefse_lda_cutoff,
        kw_cutoff       = lefse_kw_cutoff,
        wilcoxon_cutoff = lefse_wilcoxon_cutoff,
        colors          = lambda wc: GROUP_COLORS[wc.group_col]
    container:
        LEFSE_MICROBIOMEMARKER_R_CONTAINER
    resources:
        mem_mb = 6000
    script:
        f"{SCRIPTS_DIR}/run_lefse_Huttenhower.py"

rule lefse_Huttenhower_plot:
    input:
        lefse_res = TMP_DIR / "{db}/lefse_Huttenhower_results_by_{group_col}.txt"
    output:
        lda_svg       = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_LDA_by_{group_col}.svg",
        cladogram_svg = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_Cladogram_by_{group_col}.svg",
        lda_png       = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_LDA_by_{group_col}.png",
        cladogram_png = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_Cladogram_by_{group_col}.png"
    params:
        colors          = lambda wc: GROUP_COLORS[wc.group_col],
        lda_title       = "LEfSe LDA Effect Size",
        cladogram_title = "LEfSe Cladogram"
    container:
        LEFSE_PLOT_CONTAINER
    resources:
        mem_mb = 6000
    script:
        f"{SCRIPTS_DIR}/plot_lefse_res.py"
