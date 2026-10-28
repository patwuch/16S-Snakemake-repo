"""Render a `.res` file to an effect-size bar plot and a cladogram, SVG + PNG.

The bar is an LDA score for the LEfSe methods and a |log2FC| for ANCOM-BC2.
Nothing here depends on which, so the axis wording arrives as `effect_label`.

This module belongs to the plots, not to any one statistical test. Every method
whose results project onto a taxonomy tree (LEfSe Huttenhower, ANCOM-BC2)
produces a `.res` and then calls in here, so the
figures are comparable across methods instead of each library drawing in its
own idiom.

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
PLOTTERS_DIR = SCRIPT_DIR / "plotters"


def _plotter(name: str) -> str:
    p = PLOTTERS_DIR / name
    if not p.exists():
        raise FileNotFoundError(f"Plotting script not found: {p}")
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


def is_backfilled_rank(token: str, parent: str) -> bool:
    """True when write_res.R invented this rank rather than a classifier calling it.

    fill_unresolved_ranks() gives an empty rank the nearest named ancestor plus
    `_sp` (for species) or `_unresolved` (for anything else), so that every
    lineage keys its own tree node instead of piling onto a shared blank. Those
    tokens read like calls and are not: Silva138 assigns no species at all in
    some studies and still produces `s__<genus>_sp` leaves.
    """
    name = token.split("__", 1)[-1]
    parent_name = parent.split("__", 1)[-1] if parent else ""
    return bool(parent_name) and name in (f"{parent_name}_sp",
                                          f"{parent_name}_unresolved")


def deepest_significant_rank(rows) -> str:
    """Deepest taxonomic rank among nodes the method actually called significant.

    Backfilled ranks do not count. A marker set whose deepest real call is a
    genus must not report as "Species" on the strength of a synthesized leaf --
    that is a claim about the database's resolution, and it would be false.
    """
    deepest_depth = 0
    deepest_name = "Unknown"
    for fields in rows:
        if len(fields) < 3 or not fields[2]:
            continue  # no class assigned -> not significant
        ranks = fields[0].split(".")
        i = len(ranks) - 1
        while i > 0 and is_backfilled_rank(ranks[i], ranks[i - 1]):
            i -= 1
        rank_name = RANK_NAMES.get(ranks[i].split("_", 1)[0])
        if rank_name and i + 1 > deepest_depth:
            deepest_depth = i + 1
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


def render_effectsize_cladogram(res_path, barplot_svg, cladogram_svg, barplot_png, cladogram_png,
                            colors_map, db, barplot_title, cladogram_title,
                            effect_label=None, context=None, nice_ticks=False,
                            struc_zeros=None, label_boxes=False,
                            resolve_placeholders=False, legend_ncol=None) -> None:
    """Draw both plots in both formats from one `.res` file.

    `colors_map` is an ordered group -> hex mapping covering the full configured
    group universe for this axis, not just the groups that came out significant.

    `effect_label` renames the bar chart's effect-size axis for methods whose
    bars are not LDA scores (ANCOM-BC2's |log2FC|, say). Left `None` the flag is
    not passed at all, so the LEfSe methods' command lines — and therefore their
    figures — are unchanged.
   
    `context` names the subset this run was fitted on (an arm, a timepoint) and
    is prepended to the small kicker line above the title, which otherwise
    carries only the database. Without it a figure cannot say which of several
    same-titled fits it is, and that fact survives only in the report's caption.
    Left `None` the kicker is exactly what it was, so pipeline callers are
    unchanged.

    `nice_ticks` gives the bar chart a few round-numbered x ticks instead of
    LEfSe's ~10 ticks at odd values, which collide when effect sizes span a
    wide range (log2 fold changes of -8 to +9, say). Off by default, so the
    pipeline's figures are unchanged.

    `struc_zeros` is a TSV of taxa present in one group of the contrast only
    (ANCOM-BC2 structural zeros): a `#classes<TAB>first<TAB>second` line naming
    the contrast's two groups, then `lineage<TAB>group present<TAB>annotation`
    rows. The bar chart adds them below the tested bars as hatched stubs with no
    effect size; the cladogram rings their nodes. Every lineage must already be
    a node of `res_path`. Left `None`, neither plotter is told, so every
    existing figure is unchanged.

    `label_boxes` draws each keyed node's letter upright in a box of its group
    colour with a black outline, instead of as rotated text in the wedge's
    colour. Off by default.

    `resolve_placeholders` names a tip whose rank write_res.R backfilled
    ("s__COE1_sp") by the deepest rank actually assigned ("g__COE1 (no species
    call)"), in both plots. Off by default.

    `legend_ncol` sets how many columns the bar chart's legend uses. Unset, the
    plotter's default of 5 keeps every existing figure unchanged; with long
    labels a single row is wider than the chart and stretches the saved image,
    which shrinks the chart when the image is scaled to a page.
    """
    # The two plotting scripts (barplot + cladogram) re-derive the class list
    # from the .res file and sort it ALPHABETICALLY before assigning legend
    # order, color, and (for the bar plot) which side of the diverging chart
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
    # map to an underscore-only prefixed label. Both scripts split class names on
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
        # found nothing. There is no tree to draw, so take the upstream answer to
        # having nothing to plot and leave empty files behind.
        print(f"No tree nodes in {res_path}; writing empty figures.", flush=True)
        for out in (barplot_svg, cladogram_svg, barplot_png, cladogram_png):
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
    if context:
        # Subset first, then database, matching how the reports name a run in
        # prose and in every results table ("IFA+IR - Greengenes2").
        db_label = f"{context} · {db_label}"

    # A .res can be prefix-closed over a full tree and still carry no marker,
    # which is a real result rather than a failure: the tree is drawn in grey
    # and nothing is coloured. Reporting the deepest significant rank as
    # "Unknown" there reads like a lookup that fell through. Say what happened.
    if any(len(f) > 2 and f[2] for f in rows):
        footnote = (
            f"Lowest taxonomic level differentially identified: "
            f"{deepest_significant_rank(rows)}"
        )
    else:
        footnote = "No feature reached the significance cutoff"

    with tempfile.TemporaryDirectory() as _tmp:
        tmpdir = Path(_tmp)

        sz_args = []
        if struc_zeros is not None:
            sz_lines = []
            for line in Path(struc_zeros).read_text().splitlines():
                f = line.split("\t")
                if f[0] == "#classes":
                    f = [f[0]] + [sortkey[g] for g in f[1:3]]
                elif f[0]:
                    if f[1] not in sortkey:
                        raise ValueError(f"structural-zero group {f[1]!r} absent from colors_map")
                    f[1] = sortkey[f[1]]
                sz_lines.append("\t".join(f))
            (tmpdir / "struc_zeros.tsv").write_text("\n".join(sz_lines) + "\n")
            sz_args = ["--struc_zeros", str(tmpdir / "struc_zeros.tsv")]

        prefixed_res = tmpdir / "prefixed.res"
        prefixed_res.write_text("".join(
            "\t".join(sortkey.get(v, v) if i == 2 else v for i, v in enumerate(fields)) + "\n"
            for fields in rows
        ))

        # A matplotlibrc that keeps SVG text as real, searchable/selectable characters
        # instead of glyph paths. Supplied via an env var, so the (read-only /
        # encrypted) plotting scripts are untouched.
        mplrc_dir = tmpdir / "mplrc"
        mplrc_dir.mkdir()
        (mplrc_dir / "matplotlibrc").write_text("svg.fonttype: none\n")
        plot_env = {**os.environ, "MATPLOTLIBRC": str(mplrc_dir)}

        formats = [
            {
                "fmt": "svg",
                "barplot_out": barplot_svg,
                "cladogram_out": cladogram_svg,
                "label_font_size": "5",
                "barplot_extra": [],
            },
            {
                "fmt": "png",
                "barplot_out": barplot_png,
                "cladogram_out": cladogram_png,
                "label_font_size": "4",
                "barplot_extra": ["--feature_font_size", "4"],
            },
        ]

        for f in formats:
            _run([
                sys.executable, _plotter("effectsize_barplot.py"),
                str(prefixed_res),
                str(f["barplot_out"]),
                "--format", f["fmt"],
                "--dpi",    "300",
                "--left_space", "0.3",      # more space on left for labels
                *f["barplot_extra"],
                "--colors", *colors,
                "--all_feats", all_feats,
                "--class_labels", *class_labels,
                "--title", barplot_title,
                "--db_label", db_label,
                "--footnote", footnote,
                *([] if effect_label is None else ["--effect_label", effect_label]),
                *(["--nice_ticks", "1"] if nice_ticks else []),
                *sz_args,
                *(["--resolve_placeholders", "1"] if resolve_placeholders else []),
                *([] if legend_ncol is None else ["--legend_ncol", str(legend_ncol)]),
            ], env=plot_env)

            _run([
                sys.executable, _plotter("cladogram.py"),
                str(prefixed_res),
                str(f["cladogram_out"]),
                "--format", f["fmt"],
                "--dpi",    "300",
                "--colors", *colors,
                "--class_legend_font_size", "6",   # shrink group legend
                "--label_font_size",  f["label_font_size"],  # shrink feature key text
                "--labeled_start_lev", "3",        # only label from level 3 inward
                "--labeled_stop_lev",  str(depth + 1),
                "--clade_sep",        "1.5",       # more separation between clades
                "--expand_void_lev",  "1",         # expand empty levels
                "--abrv_start_lev",   "3",         # abbreviate from level 3 onward
                "--abrv_stop_lev",    str(depth),
                # A clade below this share of the circle gets its legend key on
                # the plot instead of its full name, whatever its level. Without
                # it the thin phyla each draw a ~20-character name radially at
                # the same radius and the names cross each other and every
                # abbreviation key between them.
                "--min_label_arc_deg", "20",
                "--all_feats", all_feats,
                "--class_labels", *class_labels,
                "--title", cladogram_title,
                "--db_label", db_label,
                "--footnote", footnote,
                *sz_args,
                *(["--label_boxes", "1"] if label_boxes else []),
                *(["--resolve_placeholders", "1"] if resolve_placeholders else []),
            ], env=plot_env)
