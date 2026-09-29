# ------------------------------
# Rules
# ------------------------------
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
        LEFSE_R_CONTAINER
    resources:
        mem_mb = 6000
    script:
        f"{SCRIPTS_DIR}/run_lefse_Huttenhower.py"

rule lefse_Huttenhower_plot:
    input:
        res = TMP_DIR / "{db}/lefse_Huttenhower_results_by_{group_col}.txt"
    output:
        barplot_svg   = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_LDA_by_{group_col}.svg",
        cladogram_svg = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_Cladogram_by_{group_col}.svg",
        barplot_png   = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_LDA_by_{group_col}.png",
        cladogram_png = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "lefse_Huttenhower_Cladogram_by_{group_col}.png"
    params:
        colors          = lambda wc: GROUP_COLORS[wc.group_col],
        barplot_title   = "LEfSe LDA Effect Size",
        cladogram_title = "LEfSe Cladogram"
    container:
        EFFECTSIZE_CLADOGRAM_CONTAINER
    resources:
        mem_mb = 6000
    script:
        f"{SCRIPTS_DIR}/plot_effectsize_cladogram.py"
