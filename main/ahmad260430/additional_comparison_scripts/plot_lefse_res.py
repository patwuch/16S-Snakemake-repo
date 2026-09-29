#!/usr/bin/env python3
"""Render a LEfSe .res file to LDA bar + cladogram figures (SVG and PNG).

Thin argparse wrapper around the shared scripts/effectsize_cladogram_render.py's
render_effectsize_cladogram() — that function is a plain library call with no
Snakemake coupling, so it's reused directly rather than duplicated. Must run
inside the LEfSe container (same one run_lefse_huttenhower.py used).
"""
import argparse
import sys
from pathlib import Path

REPO_SCRIPTS_DIR = Path(__file__).resolve().parents[3] / "scripts"
sys.path.insert(0, str(REPO_SCRIPTS_DIR))

from effectsize_cladogram_render import render_effectsize_cladogram  # noqa: E402


def parse_colors(pairs):
    colors = {}
    for pair in pairs:
        name, hexval = pair.split("=", 1)
        colors[name] = hexval
    return colors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--res", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--color", action="append", required=True, metavar="GROUP=HEX")
    parser.add_argument("--lda-title", default="LEfSe LDA Effect Size")
    # The renderer is method-agnostic; only the axis wording is not. Left unset
    # the flag is not forwarded, so LEfSe's own command line is unchanged.
    parser.add_argument("--effect-label", default=None)
    parser.add_argument("--cladogram-title", default="LEfSe Cladogram")
    # Which subset was fitted. Goes into the figure itself, so a panel lifted out
    # of the report still says which arm it is.
    parser.add_argument("--context", default=None)
    # Structural-zero sidecar written by structural_zero_figures.py; unset, the
    # renderer is not told and draws exactly the plain figures.
    parser.add_argument("--struc-zeros", default=None)
    # Upright, colour-boxed cladogram keys; unset, the plain layout.
    parser.add_argument("--label-boxes", action="store_true")
    # Name backfilled tips by the rank actually assigned; unset, as-is.
    parser.add_argument("--resolve-placeholders", action="store_true")
    # Bar-chart legend columns; unset, the renderer's default (one row).
    parser.add_argument("--legend-ncol", type=int, default=None)
    parser.add_argument("--lda-svg", required=True)
    parser.add_argument("--cladogram-svg", required=True)
    parser.add_argument("--lda-png", required=True)
    parser.add_argument("--cladogram-png", required=True)
    args = parser.parse_args()

    for out in (args.lda_svg, args.cladogram_svg, args.lda_png, args.cladogram_png):
        Path(out).parent.mkdir(parents=True, exist_ok=True)

    render_effectsize_cladogram(
        res_path=args.res,
        barplot_svg=args.lda_svg,
        cladogram_svg=args.cladogram_svg,
        barplot_png=args.lda_png,
        cladogram_png=args.cladogram_png,
        colors_map=parse_colors(args.color),
        db=args.db,
        barplot_title=args.lda_title,
        cladogram_title=args.cladogram_title,
        effect_label=args.effect_label,
        context=args.context,
        # Round-valued x ticks; LEfSe's default ticks collide on wide ranges.
        nice_ticks=True,
        struc_zeros=args.struc_zeros,
        label_boxes=args.label_boxes,
        resolve_placeholders=args.resolve_placeholders,
        legend_ncol=args.legend_ncol,
    )
    print(f"Wrote LEfSe figures -> {args.lda_png}, {args.cladogram_png}")


if __name__ == "__main__":
    main()
