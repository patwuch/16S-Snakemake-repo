#!/usr/bin/env python3
# Vendored from the Huttenhower LEfSe distribution (plot_res.py, shipped here
# as plugin_lefse_barplot.py until this rename) with local modifications.
# 3. ?????????
import os, sys
import matplotlib

matplotlib.use("Agg")
from pylab import *
from collections import defaultdict
import re
from lefse import *
import argparse

# colors = [
#     "#ea801c",  # D Orange
#     # "#800074",  # D Med Purple
#     "#800074",  # B Med Purple
#     "#298c8c",  # B ????????? (teal-like color
#     "#1a80bb",  # A Blue,
#     "#36b700",  # Green
#     # "#ea801c",  # A Orange
#     "#298c8c",  # ????????? (teal-like color)
#     "#f1a226",  # ????????? (orange-yellow color)
#     "#f45f74",  # Med Pink
#     "#8fd7d7",  # Light Blue
#     "#98c127",  # Med Green
#     "#ff8ca1",  # Light Pink
#     "#ffcd8e",  # Light Orange
#     "#36b700",  # Green
# ]


def read_params(args):
    parser = argparse.ArgumentParser(description="Plot results")
    parser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, help="tab delimited input file"
    )
    parser.add_argument(
        "output_file",
        metavar="OUTPUT_FILE",
        type=str,
        help="the file for the output image",
    )
    parser.add_argument(
        "--feature_font_size",
        dest="feature_font_size",
        type=int,
        default=7,
        help="the file for the output image",
    )
    parser.add_argument(
        "--format",
        dest="format",
        choices=["png", "svg", "pdf"],
        default="png",
        type=str,
        help="the format for the output file",
    )
    parser.add_argument("--dpi", dest="dpi", type=int, default=72)
    # Legend columns above the chart. The default 5 keeps every existing figure
    # unchanged; long labels in one row widen the saved canvas past the chart.
    parser.add_argument("--legend_ncol", dest="legend_ncol", type=int, default=5)
    parser.add_argument("--title", dest="title", type=str, default="")
    parser.add_argument(
        "--title_font_size", dest="title_font_size", type=str, default="12"
    )
    parser.add_argument(
        "--class_legend_font_size",
        dest="class_legend_font_size",
        type=str,
        default="10",
    )
    parser.add_argument("--width", dest="width", type=float, default=7.0)
    parser.add_argument(
        "--height",
        dest="height",
        type=float,
        default=4.0,
        help="only for vertical histograms",
    )
    parser.add_argument("--left_space", dest="ls", type=float, default=0.2)
    parser.add_argument("--right_space", dest="rs", type=float, default=0.1)
    parser.add_argument(
        "--orientation", dest="orientation", type=str, choices=["h", "v"], default="h"
    )
    parser.add_argument(
        "--autoscale", dest="autoscale", type=int, choices=[0, 1], default=1
    )
    parser.add_argument(
        "--background_color",
        dest="back_color",
        type=str,
        choices=["k", "w"],
        default="w",
        help="set the color of the background",
    )
    parser.add_argument(
        "--subclades",
        dest="n_scl",
        type=int,
        default=1,
        help="number of label levels to be dislayed (starting from the leaves, -1 means all the levels, 1 is default )",
    )
    parser.add_argument(
        "--max_feature_len",
        dest="max_feature_len",
        type=int,
        default=60,
        help="Maximum length of feature strings (def 60)",
    )
    parser.add_argument("--all_feats", dest="all_feats", type=str, default="")
    parser.add_argument(
        "--otu_only",
        dest="otu_only",
        default=False,
        action="store_true",
        help="Plot only species resolved OTUs (as opposed to all levels)",
    )
    parser.add_argument(
        "--report_features",
        dest="report_features",
        default=False,
        action="store_true",
        help="Report important features to STDOUT",
    )
    parser.add_argument(
        "--colors", dest="colors", type=str, nargs="+", help="GROUP COLORS"
    )
    parser.add_argument(
        "--footnote", dest="footnote", type=str, default="",
        help="small caption drawn at the bottom of the figure",
    )
    parser.add_argument(
        "--footnote_font_size", dest="footnote_font_size", type=str, default="6",
    )
    parser.add_argument(
        "--db_label", dest="db_label", type=str, default="",
        help="small kicker line drawn above the title (e.g. reference database name)",
    )
    parser.add_argument(
        "--nice_ticks", dest="nice_ticks", type=int, choices=[0, 1], default=0,
        help="1 = at most ~10 x ticks at round values (short labels, so no collisions). The default (0) keeps LEfSe's "
             "autoscale: ~10 ticks stepped from the axis edge at a tenth of the range, "
             "whose odd values (-6.81, -5.18, ...) collide at wide effect ranges",
    )
    parser.add_argument(
        "--effect_label", dest="effect_label", type=str,
        default="LDA SCORE (log 10)",
        help="axis label for the effect-size magnitude; defaults to LEfSe's own "
             "LDA wording so non-LEfSe methods reusing this chart can say what "
             "their bars actually measure",
    )
    parser.add_argument(
        "--resolve_placeholders", dest="resolve_placeholders", type=int, choices=[0, 1], default=0,
        help="1 = label a tip whose rank write_res.R backfilled (s__COE1_sp) by the deepest rank "
             "actually assigned (g__COE1 (no species call)). Default 0 prints the token as-is",
    )
    parser.add_argument(
        "--struc_zeros", dest="struc_zeros", type=str, default="",
        help="optional TSV of structural-zero taxa: a '#classes<TAB>first<TAB>second' line naming "
             "the contrast's two classes, then 'lineage<TAB>class where present<TAB>annotation' "
             "rows. They are present in one group only and have no effect size, so each is drawn "
             "as a hatched bar at the largest tested effect size (at an arbitrary length when "
             "nothing was tested significant) on its group's side, past that side's tested bars. "
             "Empty (default) draws exactly the plain chart",
    )
    parser.add_argument(
        "--class_labels", dest="class_labels", type=str, nargs="*", default=[],
        help="raw_class=display_text pairs; lets class values carry a sort-order-encoding "
             "prefix internally while showing clean text",
    )
    args = parser.parse_args()
    args.class_labels = dict(p.split("=", 1) for p in args.class_labels)
    args.sz, args.sz_classes = [], []
    if args.struc_zeros:
        for line in open(args.struc_zeros):
            f = line.rstrip("\n").split("\t")
            if f[0] == "#classes":
                args.sz_classes = f[1:3]
            elif f[0]:
                # Optional 4th column: "complete" (in every sample of one group) or
                # "partial" (in some). Without it every row is drawn alike.
                args.sz.append(f[:3] + [f[3] if len(f) > 3 else ""])
    return vars(args)


def read_data(input_file, output_file, otu_only, sz_classes=()):
    with open(input_file, "r") as inp:
        if not otu_only:
            rows = [
                line.strip().split()[:-1]
                for line in inp.readlines()
                if len(line.strip().split()) > 3
            ]
        else:
            rows = [
                line.strip().split()[:-1]
                for line in inp.readlines()
                if len(line.strip().split()) > 3
                and len(line.strip().split()[0].split(".")) == 8
            ]  # a feature with length 8 will have an OTU id associated with it
    classes = list(set([v[2] for v in rows if len(v) > 2]))
    if sz_classes:
        # Structural zeros are drawn for the contrast's two classes whether or not
        # any tested taxon is significant, so the class list is the contrast's.
        classes = list(sz_classes)
    if len(classes) < 1:
        print("No differentially abundant features found in " + input_file)
        # Local modification: was os.system("touch " + output_file). The output
        # path carries the study name and the grouping column, so an unquoted
        # space there made the shell create two wrong files and no declared
        # output. This branch is the common case, not a rare one.
        open(output_file, "a").close()
        sys.exit()
    data = {}
    data["rows"] = rows
    data["cls"] = classes
    return data


def plot_histo_hor(colors, path, params, data, bcl, report_features):
    cls2 = []
    if params["all_feats"] != "":
        cls2 = sorted(params["all_feats"].split(":"))
    cls = sorted(data["cls"])
    if bcl:
        data["rows"].sort(
            key=lambda ab: fabs(float(ab[3])) * (cls.index(ab[2]) * 2 - 1)
        )
    else:
        mmax = max([fabs(float(a)) for a in list(zip(*list(data["rows"])))[3]])
        data["rows"].sort(
            key=lambda ab: fabs(float(ab[3])) / mmax + (cls.index(ab[2]) + 1)
        )
    sz = params["sz"]
    if sz:
        # A structural zero is present in one group only, so it has no effect size.
        # It gets a hatched bar at the largest tested effect size (or at 1.0 when
        # nothing tested is significant), placed at the outer end of its group's
        # side, past that side's tested bars. Rows are [lineage, _, class, length,
        # "sz", annotation]; tested rows are untouched.
        n_t = len(data["rows"])
        # The first group in the colour order always takes the left side, so the
        # species and genus views of a panel line up (the plain chart picks the
        # side from whichever group's bar sorts first).
        m_sz = 1
        length = max([fabs(float(v[3])) for v in data["rows"]]) if n_t else 1.0
        side = lambda c: m_sz * (cls.index(c) * 2 - 1)
        extra = [[k, "0", c, str(length), "sz", a, s] for k, c, a, s in sorted(sz, key=lambda r: r[0])]
        # Not any(): "from pylab import *" shadows it with numpy.any, which treats a
        # generator as one truthy object.
        partial_seen = "partial" in [r[6] for r in extra]
        data["rows"] = ([r for r in extra if side(r[2]) < 0] + data["rows"]
                        + [r for r in extra if side(r[2]) > 0])
    pos = arange(len(data["rows"]))
    head = 0.75
    tail = 0.5
    ht = head + tail
    # head/tail are inch reserves, so the fractions they become shrink as the
    # figure grows with the bar count. A footnote pinned at a fixed fraction
    # therefore climbs into the x label on any tall figure. Give it inches too.
    foot = 0.45 if params["footnote"] else 0.0
    ints = max(len(pos) * 0.2, 1.5)
    height = ints + ht + foot
    fig = plt.figure(
        figsize=(params["width"], height),
        edgecolor=params["back_color"],
        facecolor=params["back_color"],
    )
    ax = fig.add_subplot(111, frame_on=False, facecolor=params["back_color"])
    ls, rs = params["ls"], 1.0 - params["rs"]
    plt.subplots_adjust(
        left=ls,
        right=rs,
        top=1 - (head * ht) / height,
        bottom=(tail * ht + foot) / height,
    )

    fig.canvas.manager.set_window_title("LDA results")

    l_align = {"horizontalalignment": "left", "verticalalignment": "baseline"}
    r_align = {"horizontalalignment": "right", "verticalalignment": "baseline"}
    added = []
    m = m_sz if sz else (1 if data["rows"][0][2] == cls[0] else -1)
    out_data = defaultdict(list)  # keep track of which OTUs result in the plot
    for i, v in enumerate(data["rows"]):
        if report_features:
            otu = (
                v[0].split(".")[7].replace("_", ".")
            )  # string replace retains format New.ReferenceOTUxxx
            score = v[3]
            otu_class = v[2]
            out_data[otu] = [score, otu_class]
        indcl = cls.index(v[2])
        lab = params["class_labels"].get(str(v[2]), str(v[2])) if str(v[2]) not in added else ""

        added.append(str(v[2]))
        col = colors[indcl % len(colors)]
        if len(cls2) > 0:
            col = colors[cls2.index(v[2]) % len(colors)]
        vv = fabs(float(v[3])) * (m * (indcl * 2 - 1)) if bcl else fabs(float(v[3]))
        if sz and len(v) > 4 and v[4] == "sz":
            if lab:
                # First bar of a group with no tested bar: key its colour solid.
                ax.barh(pos[i], 0, color=col, label=lab, height=0, edgecolor=params["fore_color"])
            if len(v) > 6 and v[6] == "partial":
                # Present in only some samples of its group: the group's own colour, a
                # sparse hatch and a dashed outline. Not a lighter colour: palettes that
                # mark irradiation by a darker shade would make a faded dark group read
                # as its light partner.
                ax.barh(pos[i], vv, align="center", height=0.8, facecolor="white", edgecolor=col,
                        hatch="//", linewidth=1.0, linestyle="--")
            else:
                ax.barh(pos[i], vv, align="center", height=0.8, facecolor="white", edgecolor=col,
                        hatch="//////", linewidth=1.0)
            ax.text(vv + (1 if vv > 0 else -1) * length / 40.0, float(i) + 0.3, v[5],
                    ha="left" if vv > 0 else "right", va="baseline",
                    size=params["feature_font_size"] * 0.85)
            continue
        ax.barh(
            pos[i],
            vv,
            align="center",
            color=col,
            label=lab,
            height=0.8,
            edgecolor=params["fore_color"],
        )
    mv = max([abs(float(v[3])) for v in data["rows"]])
    if report_features:
        print("OTU\tLDA_score\tCLass")
        for i in out_data:
            print("%s\t%s\t%s" % (i, out_data[i][0], out_data[i][1]))
    for i, r in enumerate(data["rows"]):
        indcl = cls.index(data["rows"][i][2])
        if params["n_scl"] < 0:
            rr = r[0]
        else:
            rr = ".".join(r[0].split(".")[-params["n_scl"] :])
        note = ""
        if params["resolve_placeholders"] and params["n_scl"] == 1:
            # A backfilled tip names the rank that was actually assigned, and
            # says what is missing below it, outside the italic name.
            from placeholder_names import resolve, suffix
            rr, missing = resolve(r[0])
            note = suffix(missing)
        if len(rr) > params["max_feature_len"]:
            rr = (
                rr[: params["max_feature_len"] / 2 - 2]
                + " [..]"
                + rr[-params["max_feature_len"] / 2 + 2 :]
            )

        if m * (indcl * 2 - 1) < 0 and bcl:
            match = re.search(r"(g_|s_)(.*)", rr)
            if match:
                prefix = match.group(1)  # "g_" ??? "s_"
                italic_part = match.group(2)  # ???????????????
                label_text = (
                    prefix + r"$\mathit{" + italic_part.replace("_", r"\_") + "}$"
                )
            else:
                label_text = rr
            label_text += note
            ax.text(
                mv / 40.0,
                float(i) + 0.3,
                label_text,
                l_align,
                size=params["feature_font_size"],
                color=params["fore_color"],
            )
        else:
            match = re.search(r"(g_|s_)(.*)", rr)
            if match:
                prefix = match.group(1)  # "g_" ??? "s_"
                italic_part = match.group(2)  # ???????????????
                label_text = (
                    prefix + r"$\mathit{" + italic_part.replace("_", r"\_") + "}$"
                )
            else:
                label_text = rr
            label_text += note
            ax.text(
                -mv / 40.0,
                float(i) + 0.3,
                label_text,
                r_align,
                size=params["feature_font_size"],
                color=params["fore_color"],
            )
    title_y = 1.0 + head * (1.0 - ints / (ints + ht)) * 0.8
    title_txt = ax.set_title(
        params["title"],
        size=params["title_font_size"],
        y=title_y,
        color=params["fore_color"],
    )
    db_ann = None
    if params["db_label"]:
        # Offset in points (not axes fraction) so the gap stays legible
        # regardless of figure height, which varies with feature count.
        db_ann = ax.annotate(
            params["db_label"],
            xy=(0.5, title_y),
            xycoords=("axes fraction", "axes fraction"),
            xytext=(0, float(params["title_font_size"]) * 1.6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            size=float(params["title_font_size"]) * 0.6,
            color="gray",
        )

    ax.set_yticks([])
    if sz:
        # One legend key for every hatched bar, whatever its group colour.
        if partial_seen:
            ax.barh(0, 0, height=0, facecolor="white", edgecolor="dimgray", hatch="//////",
                    label="structural zero, in every sample of one group (no log2FC)")
            ax.barh(0, 0, height=0, facecolor="white", edgecolor="dimgray", hatch="//", linestyle="--",
                    label="structural zero, in some samples of one group")
        else:
            ax.barh(0, 0, height=0, facecolor="white", edgecolor="dimgray", hatch="//////",
                    label=("structural zero (no log2FC; length = largest tested)" if n_t
                           else "structural zero (present in one group only; no log2FC)"))
    if sz and not n_t:
        # Nothing tested is on this chart: the bars mark presence only.
        ax.set_xticks([])
        ax.set_xlim(-1.35, 1.35)
        ax.set_xlabel("presence in one group only · bar length arbitrary, no log2FC")
    else:
        ax.set_xlabel(params["effect_label"])
    ax.xaxis.grid(True, color="lightgray", alpha=0.1)
    xlim = ax.get_xlim()
    if sz and not n_t:
        pass
    elif params["nice_ticks"]:
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 5, 10]))
    elif params["autoscale"]:
        ran = arange(
            0.0001, round(round((abs(xlim[0]) + abs(xlim[1])) / 10, 4) * 100, 0) / 100
        )
        if len(ran) > 1 and len(ran) < 100:
            ax.set_xticks(
                arange(
                    xlim[0],
                    xlim[1] + 0.0001,
                    min(
                        xlim[1] + 0.0001,
                        round(round((abs(xlim[0]) + abs(xlim[1])) / 10, 4) * 100, 0)
                        / 100,
                    ),
                )
            )
    ax.set_ylim((pos[0] - 1, pos[-1] + 1))
    # Row order (and hence vertical stacking) comes from `cls`/`cls2`, which are
    # both sorted so index 0 = the group earliest in the desired order. matplotlib's
    # barh plots increasing pos bottom-to-top, so without this the earliest group
    # would render at the BOTTOM — inverted relative to how the legend reads it
    # (left-to-right, earliest first). Flipping the axis (not the sort/row order,
    # which the legend and left/right bar placement still depend on) makes the
    # earliest group render at the top, so reading top-to-bottom matches the
    # desired order the same way reading the legend left-to-right does.
    ax.invert_yaxis()
    # Opt-in wrapped legend (legend_ncol other than the default 5). Three entries
    # in two columns leave one stranded, so they go in one column.
    compact = params["legend_ncol"] != 5
    ncol = params["legend_ncol"]
    if compact and ncol == 2 and len(ax.get_legend_handles_labels()[1]) == 3:
        ncol = 1
    leg = ax.legend(
        bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
        loc=3,
        ncol=ncol,
        borderaxespad=0.0,
        frameon=False,
        prop={"size": params["class_legend_font_size"]},
    )

    def get_col_attr(x):
        return hasattr(x, "set_color") and not hasattr(x, "set_facecolor")

    for o in leg.findobj(get_col_attr):
        o.set_color(params["fore_color"])
    if compact:
        # A wrapped legend is taller than the head room reserved for one row, so
        # it would run into the title on a short chart. Measure it and put the
        # title (and the kicker, which hangs off the title's anchor) above it.
        fig.canvas.draw()
        to_axes = ax.transAxes.inverted()
        top = to_axes.transform((0, leg.get_window_extent().y1))[1]
        pad = to_axes.transform((0, ax.get_window_extent().y1 + 8 / 72 * fig.dpi))[1] - 1.0
        new_y = max(title_y, top + pad)
        title_txt.set_y(new_y)
        if db_ann is not None:
            db_ann.xy = (0.5, new_y)
    for o in ax.findobj(get_col_attr):
        o.set_color(params["fore_color"])

    if params["footnote"]:
        fig.text(
            0.5,
            0.18 * foot / height,
            params["footnote"],
            ha="center",
            va="bottom",
            size=params["footnote_font_size"],
            color=params["fore_color"],
        )

    plt.savefig(
        path,
        format=params["format"],
        facecolor=params["back_color"],
        edgecolor=params["fore_color"],
        dpi=params["dpi"],
        # The title block is sized for a short chart; on a tall one its kicker
        # line lands above the canvas, as does a wrapped legend's extra rows.
        # Crop to what was drawn instead.
        **({"bbox_inches": "tight", "pad_inches": 0.2} if sz or compact else {}),
    )
    plt.close()


def plot_histo_ver(colors, path, params, data, report_features):
    cls = data["cls"]
    mmax = max([fabs(float(a)) for a in zip(*data["rows"])[1]])
    data["rows"].sort(key=lambda ab: fabs(float(ab[3])) / mmax + (cls.index(ab[2]) + 1))
    pos = arange(len(data["rows"]))
    if params["n_scl"] < 0:
        nam = [d[0] for d in data["rows"]]
    else:
        nam = [
            d[0].split(".")[-min(d[0].count("."), params["n_scl"])]
            for d in data["rows"]
        ]
    fig = plt.figure(
        edgecolor=params["back_color"],
        facecolor=params["back_color"],
        figsize=(params["width"], params["height"]),
    )
    ax = fig.add_subplot(111, facecolor=params["back_color"])
    plt.subplots_adjust(top=0.9, left=params["ls"], right=params["rs"], bottom=0.3)
    fig.canvas.manager.set_window_title("LDA results")
    l_align = {"horizontalalignment": "left", "verticalalignment": "baseline"}
    r_align = {"horizontalalignment": "right", "verticalalignment": "baseline"}
    added = []
    out_data = defaultdict(list)  # keep track of which OTUs result in the plot
    for i, v in enumerate(data["rows"]):
        if report_features:
            otu = (
                v[0].split(".")[7].replace("_", ".")
            )  # string replace retains format New.ReferenceOTUxxx
            score = v[3]
            otu_class = v[2]
            out_data[otu] = [score, otu_class]
        indcl = data["cls"].index(v[2])
        lab = str(v[2]) if str(v[2]) not in added else ""
        added.append(str(v[2]))
        col = colors[indcl % len(colors)]
        vv = fabs(float(v[3]))
        ax.bar(pos[i], vv, align="center", color=col, label=lab)
    if report_features:
        print("OTU\tLDA_score\tCLass")
        for i in out_data:
            print("%s\t%s\t%s" % (i, out_data[i][0], out_data[i][1]))
    xticks(pos, nam, rotation=-20, ha="left", size=params["feature_font_size"])
    ax.set_title(params["title"], size=params["title_font_size"])
    ax.set_ylabel(params["effect_label"])
    ax.yaxis.grid(True)
    a, b = ax.get_xlim()
    dx = float(len(pos)) / float((b - a))
    ax.set_xlim((0 - dx, max(pos) + dx))
    plt.savefig(
        path,
        format=params["format"],
        facecolor=params["back_color"],
        edgecolor=params["fore_color"],
        dpi=params["dpi"],
    )
    plt.close()


def plot_res():
    params = read_params(sys.argv)
    params["fore_color"] = "w" if params["back_color"] == "k" else "k"
    data = read_data(params["input_file"], params["output_file"], params["otu_only"], params["sz_classes"])
    if params["orientation"] == "v":
        plot_histo_ver(
            params["colors"],
            params["output_file"],
            params,
            data,
            params["report_features"],
        )
    else:
        plot_histo_hor(
            params["colors"],
            params["output_file"],
            params,
            data,
            len(data["cls"]) == 2,
            params["report_features"],
        )


if __name__ == "__main__":
    plot_res()
