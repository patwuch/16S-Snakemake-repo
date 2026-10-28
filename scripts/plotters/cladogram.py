#!/usr/bin/env python3
# Vendored from the Huttenhower LEfSe distribution (plot_cladogram.py, shipped
# here as plugin_lefse_treeplot.py until this rename) with local modifications.
# 3. ?????????
import os, sys, matplotlib, argparse, string

# from io import StringIO
matplotlib.use("Agg")
from pylab import *
from lefse import *
import numpy as np
import matplotlib.pyplot as plt
import re

# colors = [
#     "#19456e",  # Red
#     "#8f3339",  # Blue
#     # "#ea801c",  # D Orange
#     # "#800074",  # D Med Purple
#     "#800074",  # B Med Purple
#     "#298c8c",  # B ????????? (teal-like color
#     "#1a80bb",  # A Blue,
#     "#ea801c",  # A Orange
#     "#298c8c",  # ????????? (teal-like color)
#     "#f1a226",  # ????????? (orange-yellow color)
#     "#f45f74",  # Med Pink
#     "#8fd7d7",  # Light Blue
#     "#98c127",  # Med Green
#     "#ff8ca1",  # Light Pink
#     "#ffcd8e",  # Light Orange
#     "#36b700",  # Green
# ]


class CladeNode:
    def __init__(self, name, abundance, viz=True):
        self.id = name
        self.name = name.split(".")
        self.last_name = self.name[-1]
        self.abundance = abundance
        self.pos = (-1.0, -1.0)
        self.children = {}
        self.isleaf = True
        self.color = "lightgray"
        self.next_leaf = -1
        self.prev_leaf = -1
        self.viz = viz

    def __repr__(self):
        return self.last_name

    def add_child(self, node):
        self.isleaf = False
        self.children[node.__repr__()] = node

    def get_children(self):
        ck = sorted(self.children.keys())
        return [self.children[k] for k in ck]

    def get_color(self):
        return self.color

    def set_color(self, c):
        self.color = c

    def set_pos(self, pos):
        self.pos = pos


def read_params(args):
    parser = argparse.ArgumentParser(description="Cladoplot")
    parser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, help="tab delimited input file"
    )
    parser.add_argument(
        "output_file",
        metavar="OUTPUT_FILE",
        type=str,
        help="the file for the output image",
    )
    parser.add_argument("--clade_sep", dest="clade_sep", type=float, default=1.5)
    parser.add_argument("--max_lev", dest="max_lev", type=int, default=-1)
    parser.add_argument(
        "--max_point_size", dest="max_point_size", type=float, default=6.0
    )
    parser.add_argument(
        "--min_point_size", dest="min_point_size", type=float, default=1
    )
    parser.add_argument(
        "--point_edge_width", dest="markeredgewidth", type=float, default=0.25
    )
    parser.add_argument(
        "--colors", dest="colors", type=str, nargs="+", help="GROUP COLORS"
    )
    parser.add_argument(
        "--siblings_connector_width",
        dest="siblings_connector_width",
        type=float,
        default=2,
    )
    parser.add_argument(
        "--parents_connector_width",
        dest="parents_connector_width",
        type=float,
        default=0.75,
    )
    parser.add_argument(
        "--radial_start_lev", dest="radial_start_lev", type=int, default=1
    )
    parser.add_argument(
        "--labeled_start_lev", dest="labeled_start_lev", type=int, default=2
    )
    parser.add_argument(
        "--labeled_stop_lev", dest="labeled_stop_lev", type=int, default=7
    )
    parser.add_argument("--abrv_start_lev", dest="abrv_start_lev", type=int, default=3)
    # A full taxon name is drawn radially, so it sweeps across whatever sits at
    # the same radius either side of it. A clade holding a thin slice of the
    # circle cannot carry one without crossing its neighbours, however few levels
    # deep it is. Below this many degrees the name is abbreviated to its legend
    # key instead. 0 keeps the level-only behaviour.
    parser.add_argument(
        "--min_label_arc_deg", dest="min_label_arc_deg", type=float, default=0.0
    )
    parser.add_argument("--abrv_stop_lev", dest="abrv_stop_lev", type=int, default=7)
    parser.add_argument(
        "--expand_void_lev", dest="expand_void_lev", type=int, default=1
    )
    parser.add_argument(
        "--class_legend_vis", dest="class_legend_vis", type=int, default=1
    )
    parser.add_argument(
        "--colored_connector", dest="colored_connectors", type=int, default=1
    )
    parser.add_argument("--alpha", dest="alpha", type=float, default=0.2)
    parser.add_argument("--title", dest="title", type=str, default="LEfSe Cladogram")
    parser.add_argument("--sub_clade", dest="sub_clade", type=str, default="")
    parser.add_argument(
        "--title_font_size", dest="title_font_size", type=str, default="12"
    )
    parser.add_argument("--right_space_prop", dest="r_prop", type=float, default=0.1)
    parser.add_argument("--left_space_prop", dest="l_prop", type=float, default=0.1)
    parser.add_argument(
        "--label_font_size", dest="label_font_size", type=str, default="6"
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
        "--colored_labels",
        dest="col_lab",
        type=int,
        choices=[0, 1],
        default=1,
        help="draw the label with class color (1) or in black (0)",
    )
    parser.add_argument(
        "--class_legend_font_size",
        dest="class_legend_font_size",
        type=str,
        default="8",
    )
    parser.add_argument("--dpi", dest="dpi", type=int, default=72)
    parser.add_argument(
        "--format",
        dest="format",
        choices=["png", "svg", "pdf"],
        default="svg",
        type=str,
        help="the format for the output file",
    )
    parser.add_argument("--all_feats", dest="all_feats", type=str, default="")
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
        "--class_labels", dest="class_labels", type=str, nargs="*", default=[],
        help="raw_class=display_text pairs; lets class values carry a sort-order-encoding "
             "prefix internally while showing clean text",
    )
    parser.add_argument(
        "--resolve_placeholders", dest="resolve_placeholders", type=int, choices=[0, 1], default=0,
        help="1 = name a node whose rank write_res.R backfilled (s__COE1_sp) by the deepest rank "
             "actually assigned (g__COE1 (no species call)). Default 0 prints the token as-is",
    )
    parser.add_argument(
        "--label_boxes", dest="label_boxes", type=int, choices=[0, 1], default=0,
        help="1 = draw each keyed node's letter upright in a small box of its group colour with "
             "a black outline, at its usual place on the clade wedge, instead of as bare "
             "rotated text in the wedge's own colour. Default 0 keeps the plain layout",
    )
    parser.add_argument(
        "--struc_zeros", dest="struc_zeros", type=str, default="",
        help="optional TSV of structural-zero taxa ('lineage<TAB>class where present' rows; "
             "a '#classes' line is ignored). Their nodes are drawn as hollow rings in the "
             "class colour: present in one group only, so no effect size. Every lineage must "
             "be a node of the input tree. Empty (default) draws exactly the plain cladogram",
    )
    args = parser.parse_args()
    args.class_labels = dict(p.split("=", 1) for p in args.class_labels)
    args.sz = {}
    args.sz_partial = set()
    if args.struc_zeros:
        for line in open(args.struc_zeros):
            f = line.rstrip("\n").split("\t")
            if f[0] and f[0] != "#classes":
                args.sz[f[0]] = f[1]
                # Optional 4th column: "partial" rows are drawn faded.
                if len(f) > 3 and f[3] == "partial":
                    args.sz_partial.add(f[0])
    return vars(args)


def cmp_names(la, lb):
    if len(la) != len(lb):
        return False
    for p in [(a, b) for i, a in enumerate(la) for j, b in enumerate(lb) if i == j]:
        if p[0] != p[1]:
            return False
    return True


def build_tree(father, all_nodes, l, depth, viz):
    cc = [
        n
        for n in all_nodes
        if len(n.name) > len(father.name)
        and cmp_names(father.name, n.name[: len(father.name)])
    ]
    children = [n for n in cc if len(n.name) == len(father.name) + 1]
    if len(children) == 0 and l < depth - 1:  # !!!
        nc = CladeNode(father.id + "." + father.id.split(".")[-1], 1.0, viz)
        father.add_child(nc)
        children.append(nc)
    for child in children:
        build_tree(child, cc, l + 1, depth, viz)
        father.add_child(child)


def get_all_nodes(father):
    ret = [father]
    children = father.get_children()
    for c in children:
        ret += get_all_nodes(c)
    return ret


def read_data(colors, input_file, params):
    with open(input_file, "r") as inp:
        if params["sub_clade"] == "":
            rows = [
                line.strip().split()[:-1]
                for line in inp.readlines()
                if params["max_lev"] < 1
                or line.split()[0].count(".") < params["max_lev"]
            ]
        else:
            rows = [
                line.split(params["sub_clade"] + ".")[1].strip().split()[:-1]
                for line in inp.readlines()
                if (
                    params["max_lev"] < 1
                    or line.split()[0].count(".") < params["max_lev"]
                )
                and line.startswith(params["sub_clade"] + ".")
            ]
    all_names = [lin[0] for lin in rows]
    to_add = []

    abundances = [float(v) for v in list(zip(*rows))[1] if float(v) >= 0.0]
    tree = {}
    tree["classes"] = list(set([v[2] for v in rows if len(v) > 2]))
    tree["classes"].sort()
    all_nodes = [CladeNode("root." + row[0], float(row[1])) for row in rows]

    depth = max([len(n.name) for n in all_nodes])

    n2 = ["_".join(nn.name) for nn in all_nodes]
    for i, nn in enumerate(all_nodes):
        n = nn
        while "_".join(n.name[:-1]) not in n2 and len(n.name) > 1:
            n = CladeNode(".".join(n.name[:-1]), n.abundance)
            all_nodes.append(n)
            n2.append("_".join(n.name))

    cls2 = []
    if params["all_feats"] != "":
        cls2 = sorted(params["all_feats"].split(":"))
    for i, v in enumerate(rows):
        if len(v) > 2:
            if len(cls2) > 0:
                all_nodes[i].set_color(colors[cls2.index(v[2]) % len(colors)])
            else:
                if v[2].count("rgbcol") > 0:
                    ccc = [float(tt) for tt in v[2].split("_")[1:]]
                    all_nodes[i].set_color(ccc)
                else:
                    all_nodes[i].set_color(
                        colors[sorted(tree["classes"]).index(v[2]) % len(colors)]
                    )
    if params["sz"]:
        missing = set(params["sz"]) - {row[0] for row in rows}
        if missing:
            raise ValueError(f"structural-zero lineages absent from the tree: {sorted(missing)[:3]}")
        for i, v in enumerate(rows):
            c = params["sz"].get(v[0])
            if c is not None:
                if not cls2:
                    raise ValueError("--struc_zeros needs --all_feats to colour its rings")
                all_nodes[i].sz_color = colors[cls2.index(c) % len(colors)]
                all_nodes[i].sz_partial = v[0] in params["sz_partial"]
        tree["sz_classes"] = sorted(set(params["sz"].values()))
    root = CladeNode("root", -1.0)
    root.set_pos((0.0, 0.0))

    build_tree(root, all_nodes, 0, depth, params["expand_void_lev"] == 1)

    all_nodes = get_all_nodes(root)

    tree["root"] = root
    tree["max_abs"] = max(abundances)
    tree["min_abs"] = min(abundances)
    levs = []
    for i in range(depth):
        depthi = [n for n in all_nodes if len(n.name) == i + 1]
        levs.append(len(depthi))
    tree["nlev"] = levs
    return tree


def add_all_pos(father, n, distn, seps, tsep, mlev, last_leaf=-1, nc=1):
    children = father.get_children()
    leaves = True if children[0].isleaf else False
    for i, child in enumerate(children):
        if leaves:
            n += 1.0
            men = 0.5 if len(children) == 1 else 0.0
            child.set_pos(
                (
                    n * distn - men * float(distn) + tsep,
                    (len(father.name)) / float(mlev - 1),
                )
            )
            if last_leaf != -1:
                child.prev_leaf = last_leaf
                last_leaf.next_leaf = child
            last_leaf = child
        else:
            ln = n
            ltsep = tsep
            n, tsep, last_leaf = add_all_pos(
                child, n, distn, seps, tsep, mlev, last_leaf, len(children)
            )
            nn = (ln + n) * 0.5 * distn
            ssep = (ltsep + tsep) * 0.5
            if n - ln == 1:
                ssep = ltsep
            child.set_pos((nn + ssep, (len(father.name)) / float(mlev - 1)))
    tsep += seps[len(father.name) - 1]
    return n, tsep, last_leaf


def plot_points(father, params, pt_scale, ax):
    children = father.get_children()
    children.sort(key=lambda a: -int(a.get_color() == "lightgray") * a.abundance)
    x, r = father.pos[0], father.pos[1]
    for i, child in enumerate(children):
        xc, rc = plot_points(child, params, pt_scale, ax)
    if not father.viz:
        return x, r
    ps = pt_scale[0] + father.abundance / pt_scale[1] + pt_scale[0]
    col = father.get_color()
    pw = (
        params["markeredgewidth"]
        if col == "lightgray"
        else params["markeredgewidth"] * 3.0
    )
    sz_col = getattr(father, "sz_color", None)
    if sz_col is not None:
        # Structural zero: a hollow ring, never a filled marker, so it cannot be
        # read as a tested hit.
        # Partial (in some samples of its group only): the group's own colour with a
        # thin ring, not a lighter colour, which could read as another group's shade.
        partial = getattr(father, "sz_partial", False)
        ax.plot(x, r, "o", markersize=ps, markerfacecolor="white", markeredgecolor=sz_col,
                markeredgewidth=params["markeredgewidth"] * (1.5 if partial else 4.0))
        return x, r
    if x == 0 and r == 0:
        ax.plot(
            x,
            r,
            "o",
            markersize=ps,
            color=col,
            markeredgewidth=0.01,
            markeredgecolor=params["fore_color"],
        )
    else:
        ax.plot(
            x,
            r,
            "o",
            markersize=ps,
            color=col,
            markeredgewidth=pw,
            markeredgecolor=params["fore_color"],
        )

    return x, r


def plot_lines(father, params, depth, ax, xf):
    children = father.get_children()
    x, r = father.pos[0], father.pos[1]
    for i, child in enumerate(children):
        xc, rc = plot_lines(child, params, depth, ax, x)
        if i == 0:
            x_first, r_first = xc, rc
        if len(father.name) >= depth - params["radial_start_lev"]:
            col = params["fore_color"]
            lw = params["parents_connector_width"]
            if not child.viz:
                continue
            if (
                father.get_color() != "lightgray"
                and father.get_color() == child.get_color()
                and params["colored_connectors"]
            ):
                col = child.get_color()
                lw *= 2.5
            if col != params["fore_color"]:
                ax.plot([x, xc], [r, rc], "-", color=params["fore_color"], lw=lw * 1.5)
            ax.plot([x, xc], [r, rc], "-", color=col, lw=lw)

    if not father.viz or (len(children) == 1 and not children[0].viz):
        return x, r
    if len(father.name) < depth - params["radial_start_lev"]:
        col = params["fore_color"]
        lw = params["parents_connector_width"]
        if father.get_color() != "lightgray":
            f = True
            for child in children:
                if (
                    child.get_color() != father.get_color()
                    or not params["colored_connectors"]
                ):
                    f = False
                    break
            if f:
                col = father.get_color()
                lw *= 2.5
        if not (x == 0 and r == 0):
            xx = xc if len(children) > 0 else x
            if len(children) == 0:
                rc = r
            xt = x if len(children) > 1 else xx
            if col != params["fore_color"]:
                ax.plot([x, xt], [r, rc], "-", color=params["fore_color"], lw=lw * 1.5)
            ax.plot([x, xt], [r, rc], "-", color=col, lw=lw)
    if len(children) > 0 and 1 < len(father.name) < depth - params["radial_start_lev"]:
        xs = arange(x_first, xc, 0.01)
        ys = [rc for t in xs]
        ax.plot(
            xs,
            ys,
            "-",
            color=col,
            lw=params["siblings_connector_width"],
            markeredgecolor=params["fore_color"],
        )
    return x, r


def uniqueid():
    for l in string.ascii_lowercase:
        yield l
    for l in string.ascii_lowercase:
        for i in range(10):
            yield l + str(i)
    i = 0
    while True:
        yield str(i)
        i += 1


def plot_names(colors, father, params, depth, ax, u_i, seps):
    children = father.get_children()
    l = len(father.name)
    if len(children) == 0:
        if father.prev_leaf == -1 or father.next_leaf == -1:
            fr_0, fr_1 = father.pos[0], father.pos[0]
        else:
            fr_0, fr_1 = (father.pos[0] + father.prev_leaf.pos[0]) * 0.5, (
                father.pos[0] + father.next_leaf.pos[0]
            ) * 0.5
    for i, child in enumerate(children):
        fr, to = plot_names(colors, child, params, depth, ax, u_i, seps)
        if i == 0:
            fr_0 = fr
        fr_1 = to
    sz_col = getattr(father, "sz_color", None)
    if (
        father.get_color()
        != "lightgray"
        # and params["labeled_start_lev"] < l <= params["labeled_stop_lev"] + 1
    ) or sz_col is not None:
        col = father.get_color() if sz_col is None else sz_col
        group_col = col
        dd = params["labeled_stop_lev"] - params["labeled_start_lev"] + 1
        de = depth - 1
        dim = 1.0 / float(de)
        perc_ext = 0.65 if dim > 0.1 else 1.0
        clto = (de - l + 1) * dim + dim * (dd + 1 - (l - dd - 1)) * perc_ext
        clto = (de - l + 1) * dim + dim * (
            dd - (l - params["labeled_start_lev"]) + 1
        ) * perc_ext
        des = float(180.0 * (fr_0 + fr_1) / np.pi) * 0.5 - 90
        lab = ""
        txt = father.last_name
        name_tok, note = father.last_name, ""
        if params["resolve_placeholders"]:
            from placeholder_names import resolve, suffix
            real, missing = resolve(".".join(father.name))
            if missing:
                name_tok, note = real, suffix(missing)
                txt = name_tok + note
        wedge_deg = float(180.0 * (fr_1 - fr_0) / np.pi)
        too_narrow = wedge_deg < params["min_label_arc_deg"]
        if (params["abrv_start_lev"] < l or too_narrow) and l <= params["abrv_stop_lev"] + 1:
            ide = next(u_i)
            tax_level = name_tok.split("_", maxsplit=1)[0]
            tax_name = name_tok.split("_", maxsplit=1)[1].replace("_", " ")
            lab = f"{str(ide)} : {tax_level}_{tax_name}"
            txt = str(ide)
        #        ax.bar(fr_0, clto, width = fr_1-fr_0, bottom = float(l-1)/float(depth-1), alpha = params['alpha'], color=col, edgecolor=col)
        match = re.search(r"(g_|s_)(.*)", lab)
        if match:
            # ???????????????????????????????????????????????????????????????
            prefix = match.group(1)  # "g_" ??? "s_"
            italic_part = match.group(2)  # ???????????????
            label_text = (
                str(ide)
                + " : "
                + prefix
                + r"$\mathit{"
                # mathtext drops literal whitespace, and tax_name above has
                # already turned this name's underscores into spaces. Without
                # the escape, Helicobacter_C_hepaticus renders as
                # HelicobacterChepaticus.
                + italic_part.replace("_", r"\_").replace(" ", r"\ ")
                + "}$"
            )
        else:
            # ???????????????????????????????????????
            label_text = lab
        if lab:
            label_text += note

        if sz_col is not None:
            # Structural zero: the same clade wedge and key as a tested hit, but
            # hatched rather than filled, in its legend swatch too.
            ax.bar(fr_0, clto, width=fr_1 - fr_0, bottom=float(l - 1) / float(de),
                   facecolor="none", edgecolor=col, linewidth=0, alpha=0.55,
                   hatch="//" if getattr(father, "sz_partial", False) else "//////")
            partial = getattr(father, "sz_partial", False)
            # A sparse hatch can miss a small swatch entirely; the dashed outline marks
            # a partial key entry regardless.
            ax.bar(0.0, 0.0, width=0.0, bottom=0.0, facecolor="white", edgecolor=col,
                   hatch="//" if partial else "//////", linestyle="--" if partial else "-",
                   label=label_text)
        else:
            ax.bar(
                fr_0,
                clto,
                width=fr_1 - fr_0,
                bottom=float(l - 1) / float(de),
                alpha=params["alpha"],
                color=col,
                edgecolor=col,
            )
            ax.bar(
                0.0,
                0.0,
                width=0.0,
                bottom=0.0,
                alpha=1.0,
                color=col,
                edgecolor=params["fore_color"],
                # label=lab,
                label=label_text,
            )
        if l <= params["abrv_stop_lev"] + 1:
            if not params["col_lab"]:
                col = params["fore_color"]
            else:
                if col not in colors:
                    col = params["fore_color"]
                else:
                    col = colors[colors.index(col) % len(colors)]
            if params["label_boxes"] and lab:
                # Keyed node: an upright letter in a box of its group colour, so it
                # neither blends into its own (possibly hatched) wedge nor turns
                # upside down on the lower half of the circle.
                from matplotlib.colors import to_rgb
                r_, g_, b_ = to_rgb(group_col)
                ink = "black" if 0.299 * r_ + 0.587 * g_ + 0.114 * b_ > 0.55 else "white"
                drawn = ax.text(
                    (fr_0 + fr_1) * 0.5,
                    clto + float(l - 1) / float(de) - dim * perc_ext / 2.0,
                    txt,
                    size=params["label_font_size"],
                    rotation=0,
                    ha="center",
                    va="center",
                    color=ink,
                    fontweight="bold",
                    zorder=6,
                    bbox=dict(boxstyle="square,pad=0.25", facecolor=group_col, edgecolor="black",
                              linewidth=0.5),
                )
            else:
                drawn = ax.text(
                    (fr_0 + fr_1) * 0.5,
                    clto + float(l - 1) / float(de) - dim * perc_ext / 2.0,
                    txt,
                    size=params["label_font_size"],
                    rotation=des,
                    ha="center",
                    va="center",
                    color=col,
                )
            params["_labels"].append(drawn)
    return fr_0, fr_1


def declutter_labels(fig, labels, max_passes=60, step=0.015, max_shift=0.6):
    """Slide overlapping labels outward along the radius until they clear.

    Sibling clades holding a thin slice of the circle put their labels at the
    same radius a couple of degrees apart, so the keys land on top of each other
    and none of them can be read. By the time the labels exist the angles are
    fixed by the tree, so radius is the only axis left to spend.

    Each pass moves the outer label of every overlapping pair one step further
    out and re-measures. `max_shift` caps the total displacement so a label in a
    hopelessly crowded sector settles for being close rather than drifting off
    the figure.
    """
    if not labels:
        return
    home = [t.get_position()[1] for t in labels]
    for _ in range(max_passes):
        fig.canvas.draw()
        boxes = [t.get_window_extent() for t in labels]
        moved = False
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if not boxes[i].overlaps(boxes[j]):
                    continue
                out = i if labels[i].get_position()[1] >= labels[j].get_position()[1] else j
                theta, r = labels[out].get_position()
                if r + step - home[out] > max_shift:
                    continue
                labels[out].set_position((theta, r + step))
                moved = True
        if not moved:
            return


def draw_tree(colors, out_file, tree, params):
    plt_size = 10
    nlev = tree["nlev"]
    pt_scale = (
        params["min_point_size"],
        max(
            1.0,
            ((tree["max_abs"] - tree["min_abs"]))
            / (params["max_point_size"] - params["min_point_size"]),
        ),
    )
    depth = len(nlev)
    sep = (2.0 * np.pi) / float(nlev[-1])
    seps = [
        params["clade_sep"] * sep / float(depth - i + 1)
        for i in range(1, len(tree["nlev"]) + 1)
    ]
    totseps = sum([s * nlev[i] for i, s in enumerate(seps[:-1])])
    clade_sep_err = True if totseps > np.pi else False
    while totseps > np.pi:
        params["clade_sep"] *= 0.75
        seps = [
            params["clade_sep"] * sep / (float(depth - i + 1) * 0.25)
            for i in range(1, len(tree["nlev"]) + 1)
        ]
        totseps = sum([s * nlev[i] for i, s in enumerate(seps[:-1])])
    if clade_sep_err:
        print("clade_sep parameter too large, lowered to", params["clade_sep"])

    fig = plt.figure(edgecolor=params["back_color"], facecolor=params["back_color"])
    ax = fig.add_subplot(
        111, polar=True, frame_on=False, facecolor=params["back_color"]
    )
    # No margins are reserved for the legends. They are placed outside the drawn
    # cladogram once it has been measured, and the tight crop at save grows the
    # canvas to hold them, so --left_space_prop / --right_space_prop are unused.
    ax.grid(False)
    xticks([])
    yticks([])

    ds = (2.0 * np.pi - totseps) / float(nlev[-1])

    add_all_pos(tree["root"], 0.0, ds, seps, 0.0, depth)

    plot_lines(tree["root"], params, depth, ax, 0)
    plot_points(tree["root"], params, pt_scale, ax)
    params["_labels"] = []
    plot_names(colors, tree["root"], params, depth, ax, uniqueid(), seps)
    declutter_labels(fig, params["_labels"])
    r = np.arange(0, 3.0, 0.01)
    theta = 2 * np.pi * r

    def get_col_attr(x):
        return hasattr(x, "set_color") and not hasattr(x, "set_facecolor")

    ax.set_title(
        params["title"], size=params["title_font_size"], color=params["fore_color"]
    )
    db_kicker = None
    if params["db_label"]:
        db_kicker = ax.annotate(
            params["db_label"],
            xy=(0.5, 1.0),
            xycoords=("axes fraction", "axes fraction"),
            xytext=(0, float(params["title_font_size"]) * 2.2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            size=float(params["title_font_size"]) * 0.6,
            color="gray",
        )

    # Everything drawn so far -- tree, outside clade labels, title, kicker -- is
    # the cladogram. Measure it once, before any legend exists, and put the
    # legends and footnote outside this box. Nothing is ever pulled back inside
    # the canvas: the tight crop at save grows the image to hold them instead,
    # so no key size can land on the drawing.
    fig.canvas.draw()
    clad = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(
        fig.transFigure.inverted()
    )
    fw, fh = fig.get_size_inches()
    pad_x, pad_y = 0.1 / fw, 0.1 / fh

    h, l = ax.get_legend_handles_labels()
    taxa_leg = None
    if len(l) > 0:
        # Each column allows at most 35 species (rows)
        ncol = len(l) // 35 + 1
        taxa_leg = ax.legend(
            bbox_to_anchor=(clad.x1 + pad_x, clad.y1),
            bbox_transform=fig.transFigure,
            frameon=False,
            loc="upper left",
            borderaxespad=0.0,
            prop={"size": params["label_font_size"]},
            ncol=ncol,
        )
        if taxa_leg != None:
            gca().add_artist(taxa_leg)
            for o in taxa_leg.findobj(get_col_attr):
                o.set_color(params["fore_color"])

    cll = (
        sorted(tree["classes"])
        if params["all_feats"] == ""
        else sorted(params["all_feats"].split(":"))
    )
    nll = [
        ax.bar(0.0, 0.0, width=0.0, bottom=0.0, color=colors[i % len(colors)], label=c)
        for i, c in enumerate(cll)
        if c in tree["classes"]
    ]
    cl = [params["class_labels"].get(c, c) for c in cll if c in tree["classes"]]
    if params["sz"]:
        # Key every group that has a ring, even one with no tested hit, then the
        # ring itself.
        shown = set(tree["classes"]) | set(tree["sz_classes"])
        nll = [ax.bar(0.0, 0.0, width=0.0, bottom=0.0, color=colors[i % len(colors)], label=c)
               for i, c in enumerate(cll) if c in shown]
        cl = [params["class_labels"].get(c, c) for c in cll if c in shown]
        from matplotlib.lines import Line2D
        nll.append(Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor="white",
                          markeredgecolor="dimgray", markeredgewidth=1.5))
        if params["sz_partial"]:
            cl.append("hollow ring / hatched: structural zero, in every sample of one group")
            nll.append(Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor="white",
                              markeredgecolor="dimgray", markeredgewidth=0.6))
            cl.append("thin ring / sparse hatch: structural zero, in some samples of one group")
        else:
            cl.append("hollow ring / hatched: structural zero (no log2FC)")

    class_leg = None
    if params["class_legend_vis"]:
        # Left of the cladogram, level with its bottom edge.
        class_leg = l2 = legend(
            nll,
            cl,
            loc="lower right",
            bbox_to_anchor=(clad.x0 - pad_x, clad.y0),
            bbox_transform=fig.transFigure,
            prop={"size": params["class_legend_font_size"]},
            frameon=False,
        )
        if l2 != None:
            for o in l2.findobj(get_col_attr):
                o.set_color(params["fore_color"])
    footnote = None
    if params["footnote"]:
        # Under the cladogram, below its lowest outside label.
        footnote = fig.text(
            (clad.x0 + clad.x1) / 2,
            clad.y0 - pad_y,
            params["footnote"],
            ha="center",
            va="top",
            size=params["footnote_font_size"],
            color=params["fore_color"],
        )

    # The legends and footnote may lie outside the canvas; crop to what was
    # actually drawn, which grows the image to hold them. Both legends have to be named here: a legend parented to the axes is not
    # part of the figure's own tight bbox, so cropping without them slices the
    # taxa key off mid-word. The footnote goes in the same list because
    # bbox_extra_artists REPLACES matplotlib's default artist set rather than
    # extending it, and figure-level text is in that default. A run with no
    # markers draws no taxa key, so nothing else holds the crop open and the
    # footnote was being cut mid-word.
    # The kicker is in the same position for a subtler reason: Figure.get_tightbbox
    # forwards bbox_extra_artists down to each Axes, so an axes that is handed a
    # list stops considering its own children. Only the title is special-cased
    # back in, which is why the title survived the crop and the kicker above it
    # did not.
    extra = [a for a in (taxa_leg, class_leg, footnote, db_kicker) if a is not None]

    plt.savefig(
        out_file,
        format=params["format"],
        facecolor=params["back_color"],
        edgecolor=params["fore_color"],
        dpi=params["dpi"],
        bbox_inches="tight",
        bbox_extra_artists=extra,
    )
    plt.close()


def plot_cladogram():
    params = read_params(sys.argv)
    params["fore_color"] = "w" if params["back_color"] == "k" else "k"
    clad_tree = read_data(params["colors"], params["input_file"], params)
    draw_tree(params["colors"], params["output_file"], clad_tree, params)


if __name__ == "__main__":

    plot_cladogram()
