import os
# Send this job's stdout/stderr, subprocesses included, to its Snakemake log.
_log = os.open(snakemake.log[0], os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
os.dup2(_log, 1)
os.dup2(_log, 2)

from effectsize_cladogram_render import render_effectsize_cladogram

render_effectsize_cladogram(
    res_path        = snakemake.input.res,
    barplot_svg     = snakemake.output.barplot_svg,
    cladogram_svg   = snakemake.output.cladogram_svg,
    barplot_png     = snakemake.output.barplot_png,
    cladogram_png   = snakemake.output.cladogram_png,
    colors_map      = snakemake.params.colors,
    db              = snakemake.wildcards.db,
    barplot_title   = snakemake.params.barplot_title,
    cladogram_title = snakemake.params.cladogram_title,
    # Optional: only methods whose bars are not LDA scores set it, and leaving
    # it None keeps the LEfSe barplot command line byte-identical.
    effect_label    = getattr(snakemake.params, "effect_label", None),
)
