from lefse_render import render_lefse_figures

render_lefse_figures(
    res_path       = snakemake.input.lefse_res,
    lda_svg        = snakemake.output.lda_svg,
    cladogram_svg  = snakemake.output.cladogram_svg,
    lda_png        = snakemake.output.lda_png,
    cladogram_png  = snakemake.output.cladogram_png,
    colors_map     = snakemake.params.colors,
    db             = snakemake.wildcards.db,
    lda_title       = snakemake.params.lda_title,
    cladogram_title = snakemake.params.cladogram_title,
)
