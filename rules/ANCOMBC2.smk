rule ancombc2_run:
    input:
        feature_table = DA_FEATURE_TABLE,
        metadata = STUDY_DIR / "metadata.tsv",
        taxonomy = DA_TAXONOMY
    output:
        output_file = TMP_DIR / "{db}/ANCOMBC2_results_by_{group_col}.tsv",
        # The shared .res interface, so ANCOM-BC2's FDR-controlled hits can be
        # drawn as an effect-size bar chart plus cladogram and compared against
        # LEfSe's permissive ones. The .tsv beside it is the raw ancombc2()
        # table, kept for inspection; nothing downstream reads it.
        res_file    = TMP_DIR / "{db}/ANCOMBC2_results_by_{group_col}.res"
    container:
        QIIME_CONTAINER
    params:
        group_col            = "{group_col}",
        # Picks the per-database collapse level, same policy LEfSe applies.
        db                   = "{db}",
        min_sample_presence  = ancombc2_min_sample_presence,
        prv_cut              = ancombc2_prv_cut,
        lib_cut              = ancombc2_lib_cut,
        p_adj_method         = ancombc2_p_adj_method,
        alpha                = ancombc2_alpha,
        pseudo_sens          = ancombc2_pseudo_sens,
        struc_zero           = ancombc2_struc_zero,
        neg_lb               = ancombc2_neg_lb,
        precollapsed         = lambda wc: region == "region_V3V4" and wc.db == "Silva138"
    script:
        f"{SCRIPTS_DIR}/run_ancombc2.R"

# Drawn through the renderer the LEfSe methods also use. The bars are log2 fold
# changes, not LDA scores, which is why the file name does not say LDA.
rule ancombc2_effectsize_cladogram_plot:
    input:
        res = TMP_DIR / "{db}/ANCOMBC2_results_by_{group_col}.res"
    output:
        barplot_svg   = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "ANCOMBC2_EffectSize_by_{group_col}.svg",
        cladogram_svg = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "ANCOMBC2_Cladogram_by_{group_col}.svg",
        barplot_png   = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "ANCOMBC2_EffectSize_by_{group_col}.png",
        cladogram_png = DIFFERENTIAL_ABUNDANCE_DIR / "{db}" / "ANCOMBC2_Cladogram_by_{group_col}.png"
    params:
        colors          = lambda wc: GROUP_COLORS[wc.group_col],
        barplot_title   = "ANCOM-BC2 Effect Size",
        cladogram_title = "ANCOM-BC2 Cladogram",
        effect_label    = "|log2FC| (mdFDR-controlled)"
    container:
        EFFECTSIZE_CLADOGRAM_CONTAINER
    resources:
        mem_mb = 4000
    script:
        f"{SCRIPTS_DIR}/plot_effectsize_cladogram.py"