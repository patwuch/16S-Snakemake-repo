"""Render a LEfSe `.res` file to an LDA bar plot and a cladogram, SVG + PNG.

Every LEfSe method in this pipeline (Huttenhower, MicrobiomeMarker) produces a
`.res` and then calls in here, so the figures are comparable across methods
instead of each library drawing in its own idiom.

The `.res` file is the whole interface: a prefix-closed table over the taxonomy
tree, one TAB-separated row per node, written by whichever method produced it.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

# Only used for the small kicker line drawn above the plot title. An unlisted
# database is still a perfectly runnable analysis, so fall back to its raw name
# rather than failing the render.
DB_DISPLAY = {
    "Silva138":    "SILVA 138",
    "Greengenes2": "Greengenes2",
}

RANK_NAMES = {
    "d": "Domain", "p": "Phylum", "c": "Class",
    "o": "Order", "f": "Family", "g": "Genus", "s": "Species",
}

SCRIPT_DIR = Path(__file__).resolve().parent
LEFSE_DIR = SCRIPT_DIR / "lefse"


def _lefse(name: str) -> str:
    p = LEFSE_DIR / name
    if not p.exists():
        raise FileNotFoundError(f"LEfSe script not found: {p}")
    return str(p)


def _run(cmd: list, env=None) -> None:
    print("Running:", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True, env=env)


def _read_rows(res_path):
    with open(res_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                yield line.split("\t")


def deepest_significant_rank(rows) -> str:
    """Deepest taxonomic rank among nodes the method actually called significant."""
    deepest_depth = 0
    deepest_name = "Unknown"
    for fields in rows:
        if len(fields) < 3 or not fields[2]:
            continue  # no class assigned -> not significant
        ranks = fields[0].split(".")
        rank_name = RANK_NAMES.get(ranks[-1].split("_", 1)[0])
        if rank_name and len(ranks) > deepest_depth:
            deepest_depth = len(ranks)
            deepest_name = rank_name
    return deepest_name


def max_node_depth(rows) -> int:
    """Depth of the deepest node in the tree, in dot-separated rank units.

    The cladogram's labelling/abbreviation cutoffs used to be a per-database
    constant, on the theory that Greengenes2 resolves species and SILVA 138 does
    not. That stopped being true once the R methods joined: they read the
    genus-collapsed BIOM and are genus-capped whatever the database. The tree
    that was actually written is the only honest source for its own depth.
    """
    return max((len(f[0].split(".")) for f in rows), default=1)


def render_lefse_figures(res_path, lda_svg, cladogram_svg, lda_png, cladogram_png,
                         colors_map, db, lda_title, cladogram_title) -> None:
    """Draw both plots in both formats from one `.res` file.

    `colors_map` is an ordered group -> hex mapping covering the full configured
    group universe for this axis, not just the groups that came out significant.
    """
    # The LEfSe plotting scripts (barplot + cladogram) re-derive the class list
    # from the .res file and sort it ALPHABETICALLY before assigning legend
    # order, color, and (for the LDA barplot) which side of the diverging chart
    # each group's bars land on. Row order in the .res file is ignored, so
    # reordering rows does nothing. To control that ordering we make the
    # alphabetical order EQUAL the desired order by prefixing each class label
    # with a sort key ("zzsort00_", "zzsort01_", ...). The prefix goes on a temp
    # copy only: the persisted .res is a real pipeline output and keeps plain
    # group names. The TEXT the scripts draw comes from `--class_labels`, a
    # raw->display lookup, so the figure never shows the prefix, for either SVG
    # or PNG. (PNG's rasterized text can't be edited after the fact the way SVG
    # text can, which is why display text has to be supplied at render time
    # rather than stripped out afterward.)
    order = list(colors_map.keys())              # desired legend order
    colors = [colors_map[g] for g in order]      # colors in DESIRED order
    #   (because the prefixed labels now sort into `order`, the plotting scripts
    #    assign colors[i] to order[i] — so pass colors in desired order, not sorted.)

    # Both the original (space) and QIIME 2 (underscore) forms of each group name
    # map to an underscore-only prefixed label. LEfSe scripts split class names on
    # whitespace internally, so spaces in class labels cause parse errors;
    # underscores are safe.
    sortkey = {}
    for i, g in enumerate(order):
        g_norm = g.replace(" ", "_")
        label = f"zzsort{i:02d}_{g_norm}"
        sortkey[g] = label
        sortkey[g_norm] = label

    all_feats = ":".join(sortkey[g] for g in order)
    class_labels = [f"{sortkey[g]}={g}" for g in order]

    rows = list(_read_rows(res_path))
    if not rows:
        # A .res with no nodes means the comparison was degenerate, not that it
        # found nothing. There is no tree to draw, so take LEfSe's own answer to
        # having nothing to plot and leave empty files behind.
        print(f"No tree nodes in {res_path}; writing empty figures.", flush=True)
        for out in (lda_svg, cladogram_svg, lda_png, cladogram_png):
            Path(out).touch()
        return

    unknown = sorted({f[2] for f in rows
                      if len(f) > 2 and f[2] and f[2] not in sortkey})
    if unknown:
        raise ValueError(
            f"Enrichment class(es) in {res_path} absent from colors_map: {unknown}.\n"
            f"colors_map keys: {sorted(colors_map)}"
        )

    depth = max_node_depth(rows)
    db_label = DB_DISPLAY.get(db, db)
    footnote = (
        f"Lowest taxonomic level differentially identified: "
        f"{deepest_significant_rank(rows)}"
    )

    with tempfile.TemporaryDirectory() as _tmp:
        tmpdir = Path(_tmp)

        prefixed_res = tmpdir / "prefixed.res"
        prefixed_res.write_text("".join(
            "\t".join(sortkey.get(v, v) if i == 2 else v for i, v in enumerate(fields)) + "\n"
            for fields in rows
        ))

        # A matplotlibrc that keeps SVG text as real, searchable/selectable characters
        # instead of glyph paths. Supplied via an env var, so the (read-only /
        # encrypted) lefse scripts are untouched.
        mplrc_dir = tmpdir / "mplrc"
        mplrc_dir.mkdir()
        (mplrc_dir / "matplotlibrc").write_text("svg.fonttype: none\n")
        plot_env = {**os.environ, "MATPLOTLIBRC": str(mplrc_dir)}

        formats = [
            {
                "fmt": "svg",
                "lda_out": lda_svg,
                "cladogram_out": cladogram_svg,
                "label_font_size": "5",
                "barplot_extra": [],
            },
            {
                "fmt": "png",
                "lda_out": lda_png,
                "cladogram_out": cladogram_png,
                "label_font_size": "4",
                "barplot_extra": ["--feature_font_size", "4"],
            },
        ]

        for f in formats:
            _run([
                sys.executable, _lefse("plugin_lefse_barplot.py"),
                str(prefixed_res),
                str(f["lda_out"]),
                "--format", f["fmt"],
                "--dpi",    "300",
                "--left_space", "0.3",      # more space on left for labels
                *f["barplot_extra"],
                "--colors", *colors,
                "--all_feats", all_feats,
                "--class_labels", *class_labels,
                "--title", lda_title,
                "--db_label", db_label,
                "--footnote", footnote,
            ], env=plot_env)

            _run([
                sys.executable, _lefse("plugin_lefse_treeplot.py"),
                str(prefixed_res),
                str(f["cladogram_out"]),
                "--format", f["fmt"],
                "--dpi",    "300",
                "--colors", *colors,
                "--left_space_prop",  "0.15",
                "--right_space_prop", "0.45",
                "--class_legend_font_size", "6",   # shrink group legend
                "--label_font_size",  f["label_font_size"],  # shrink feature key text
                "--labeled_start_lev", "3",        # only label from level 3 inward
                "--labeled_stop_lev",  str(depth + 1),
                "--clade_sep",        "1.5",       # more separation between clades
                "--expand_void_lev",  "1",         # expand empty levels
                "--abrv_start_lev",   "3",         # abbreviate from level 3 onward
                "--abrv_stop_lev",    str(depth),
                "--all_feats", all_feats,
                "--class_labels", *class_labels,
                "--title", cladogram_title,
                "--db_label", db_label,
                "--footnote", footnote,
            ], env=plot_env)
