# DB-agnostic custom taxa barplot plotting. Wildcarded on {db}, so it must be
# included unconditionally rather than living inside a per-DB rules file — only
# one of Greengenes2.smk/Silva138.smk is included depending on REFERENCE_DB, but
# every configured DB still needs this rule.
rule plot_taxa_barplot:
    input:
        table_biom = TABLES_DIR / "study-seqs.biom",
        taxonomy_tsv = TABLES_DIR / "exported-taxonomy" / "{db}_taxonomy.tsv"
    output:
        plot_samples     = str(TAXA_BARPLOT_DIR / "{db}" / "taxa_barplot_{taxa_level}_by_{factor}_samples.svg"),
        plot_groups      = str(TAXA_BARPLOT_DIR / "{db}" / "taxa_barplot_{taxa_level}_by_{factor}_groups.svg"),
        plot_samples_png = str(TAXA_BARPLOT_DIR / "{db}" / "taxa_barplot_{taxa_level}_by_{factor}_samples.png"),
        plot_groups_png  = str(TAXA_BARPLOT_DIR / "{db}" / "taxa_barplot_{taxa_level}_by_{factor}_groups.png")
    message:
        "Plotting taxa barplot for level={wildcards.taxa_level}, factor={wildcards.factor}, db={wildcards.db}"
    params:
        group_by     = "{factor}",
        group_order  = lambda wc: GROUP_ORDERS[wc.factor],
        taxa_level   = "{taxa_level}",
        database     = "{db}",
        top_n_taxa_shown_on_barplot = top_n_taxa,
        metadata_tsv = metadata_path,
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    log:
        STUDY_DIR / "logs/plot_taxa_barplot/{db}_{factor}_{taxa_level}.log"
    script:
        SCRIPTS_DIR / "plot_taxa_barplot_levels.py"
