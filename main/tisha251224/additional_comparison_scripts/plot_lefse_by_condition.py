#!/usr/bin/env python3
"""Re-draw the pipeline's LEfSe figures, whole and split by Condition, without refitting.

LEfSe was run once across all six Group levels; each marker is enriched in one
group against the other five. This draws that result as the pipeline does, then
twice more: once with only the markers enriched in an IR group, once with only
those enriched in a Control group. Every other node keeps its place in the tree,
uncoloured.

Group colours and order come from metadata.tsv (PrimaryColor, PrimaryOrder).

Also writes tables/<db>_lefse_Huttenhower_markers_by_Group.tsv, as the pipeline
does, since this study's LEfSe ran before that table existed.

Must run inside the LEfSe container, like the pipeline's lefse_Huttenhower_plot.
"""
import csv
import sys
from pathlib import Path

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY.parents[1] / "scripts"))

from effectsize_cladogram_render import render_effectsize_cladogram  # noqa: E402
from lefse_marker_table import write_lefse_marker_table  # noqa: E402

with open(STUDY / "metadata.tsv") as fh:
    groups = {r["Group"]: (int(r["PrimaryOrder"]), r["PrimaryColor"])
              for r in csv.DictReader(fh, delimiter="\t")}
COLORS = {g: c for g, (_, c) in sorted(groups.items(), key=lambda kv: kv[1][0])}
SUBSETS = {"IR": "Irradiated", "Control": "Unirradiated"}

for db in ("Greengenes2", "Silva138"):
    res = STUDY / "tmp" / db / "lefse_Huttenhower_results_by_Group.txt"
    out = STUDY / "plots" / "differential_abundance" / db
    # The table the pipeline's lefse_Huttenhower_table rule writes.
    write_lefse_marker_table(res, STUDY / "tables" / f"{db}_lefse_Huttenhower_markers_by_Group.tsv",
                             list(COLORS))

    # The pipeline's own figures, with the pipeline's own options.
    render_effectsize_cladogram(
        res_path        = res,
        barplot_svg     = out / "lefse_Huttenhower_LDA_by_Group.svg",
        cladogram_svg   = out / "lefse_Huttenhower_Cladogram_by_Group.svg",
        barplot_png     = out / "lefse_Huttenhower_LDA_by_Group.png",
        cladogram_png   = out / "lefse_Huttenhower_Cladogram_by_Group.png",
        colors_map      = COLORS,
        db              = db,
        barplot_title   = "LEfSe LDA Effect Size",
        cladogram_title = "LEfSe Cladogram",
    )

    for cond, label in SUBSETS.items():
        colors = {g: c for g, c in COLORS.items() if g.endswith(f"_{cond}")}
        sub_res = res.with_name(f"{res.stem}_{cond}.txt")
        with open(res) as src, open(sub_res, "w") as dst:
            for line in src:
                f = line.rstrip("\n").split("\t")
                if len(f) > 2 and f[2] and f[2] not in colors:
                    f[2:] = ["", "", "-"]  # LEfSe's own non-marker row
                dst.write("\t".join(f) + "\n")
        render_effectsize_cladogram(
            res_path        = sub_res,
            barplot_svg     = out / f"lefse_Huttenhower_LDA_by_Group_{cond}.svg",
            cladogram_svg   = out / f"lefse_Huttenhower_Cladogram_by_Group_{cond}.svg",
            barplot_png     = out / f"lefse_Huttenhower_LDA_by_Group_{cond}.png",
            cladogram_png   = out / f"lefse_Huttenhower_Cladogram_by_Group_{cond}.png",
            colors_map      = colors,
            db              = db,
            barplot_title   = "LEfSe LDA Effect Size",
            cladogram_title = "LEfSe Cladogram",
            context         = label,
            # One row of three long labels is wider than the chart and clips.
            legend_ncol     = 1,
        )
