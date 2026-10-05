#!/usr/bin/env python3
"""HTML for report.html's MaAsLin3 table body and the whole-design figures, from the tables.

Prints two fragments separated by a line "<!-- panels -->": the <tbody> of the
DA table, and the "Both together" figures (clustered factorial heatmap and
faceted forest plot from plot_factorial_views.py, species and genus views, with
captions naming the robust calls). The prose around them is written by hand;
check_report.py verifies both.

Usage, inside the qiime2 container:
    python3 gen_da_html.py > fragments.html
"""
import html as H
import os
import re

import pandas as pd

from da_rows import Q_INT, rows, species_named

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
TAB = os.path.join(OUT, "tables")

meta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t", index_col=0, keep_default_na=False)
COLOR = dict(zip(meta.Group, meta.PrimaryColor))
MINUS = "−"
LONG = {"Cellulose": "cellulose", "Inulin": "inulin", "RS": "resistant starch"}


def num(v, plus=True):
    return (("+" if plus else "") + f"{v:.2f}") if v >= 0 else f"{MINUS}{-v:.2f}"


def taxon_html(t):
    """Deepest assigned rank, italicised the way the report names taxa."""
    parts = [p for p in t.split(";") if not p.endswith("__")]
    rank, name = parts[-1][0], parts[-1][3:]
    clean = lambda s: re.sub(r"_\d+$", "", s)
    if rank == "s":
        genus, _, epithet = name.partition(" ")
        gm = re.match(r"(.+?)(_[A-Z])?(?:_\d+)?$", genus)
        g, gsuf = gm.group(1), gm.group(2) or ""
        epithet = clean(epithet)
        if re.match(r"sp\d+", epithet):
            return f'<em class="taxon">{g}</em>{gsuf} {epithet}'
        em = re.match(r"(.+?)(_[A-Z])?$", epithet)
        if gsuf:
            return f'<em class="taxon">{g}</em>{gsuf} <em class="taxon">{em.group(1)}</em>{em.group(2) or ""}'
        return f'<em class="taxon">{g} {em.group(1)}</em>{em.group(2) or ""}'
    word = {"g": "genus", "f": "family", "o": "order"}[rank]
    below = {"g": "species", "f": "genus", "o": "family"}[rank]
    return f'{word} <em class="taxon">{clean(name)}</em> (no {below} call)'


def cell(v):
    return f'<td class="n">{v}</td>'


def table_body(rs):
    out, q = [], None
    for r in rs:
        head = ""
        if r["question"] != q:
            q = r["question"]
            k = sum(x["question"] == q for x in rs)
            head = f'<tr class="grp"><td rowspan="{k}" class="term q">{H.escape(q, quote=False)}</td>'
        else:
            head = "<tr>"
        sw = "" if r["question"] == Q_INT else f'<span class="swatch" style="background:{COLOR[r["first"]]}"></span>'
        ud = "—" if not r["n"] else f'{"+" + format(r["max_up"], ".2f") if r["max_up"] is not None else "—"} / ' \
                                    f'{MINUS + format(r["max_down"], ".2f") if r["max_down"] is not None else "—"}'
        mq = "—" if r["min_q"] is None else f'{r["min_q"]:.1e}'
        chip, kind = r["chip"]
        calls = f'{r["n"]} ({r["n_robust"]})'
        out.append(f'{head}<td>{sw}{r["label"]}</td>{cell(r["fitted"])}{cell(calls)}'
                   f'{cell(ud)}{cell(mq)}{cell(num(r["median"]))}<td><span class="chip {kind}">{chip}</span></td></tr>')
    return "<tbody>\n" + "\n".join(out) + "\n</tbody>"


PRE = "MaAsLin3_factorial_Greengenes2"


def call_phrase(c):
    return f'{taxon_html(c.taxon)} ({num(c.coef)}, q = {c.q:.2g})'


def robust_list(rs, keep=lambda t: True):
    """Robust calls, grouped by row: "Inulin vs Cellulose: A (+x, q = y); B (...)"."""
    parts = []
    for r in rs:
        rob = r["calls"][r["calls"].robust]
        rob = rob[[keep(t) for t in rob.taxon]]
        if len(rob):
            parts.append(f'{r["label"]}{"" if r["question"] != Q_INT else " (interaction)"}: '
                         + "; ".join(call_phrase(c) for c in rob.itertuples()))
    return parts


def panel(i, title, stem, alt, cap_sp, cap_gn):
    gn = stem.replace(PRE, f"{PRE}_gplus")
    return f'''<figure class="panel">
  <div class="panel-head"><span class="panel-title">{title}</span>
    <span class="seg"><input type="radio" name="v{i}" id="v{i}s" value="sp" checked><label for="v{i}s">Species</label><input type="radio" name="v{i}" id="v{i}g" value="gn"><label for="v{i}g">Genus and above</label></span></div>
    <div class="view sp">
      <div class="one">
      <img src="figures/{stem}.png" alt="{alt}">
      </div>
      <figcaption>{cap_sp}</figcaption>
    </div>
    <div class="view gn">
      <div class="one">
      <img src="figures/{gn}.png" alt="{alt}, taxa named only to genus or family">
      </div>
      <figcaption>{cap_gn}</figcaption>
    </div>
</figure>'''


def factorial_views(rs):
    taxa = sorted({t for r in rs for t in r["calls"].taxon})
    genus = [t for t in taxa if not species_named(t)]
    at6 = all(r["calls"].call_at_6.all() for r in rs if r["n"])
    rob = robust_list(rs)
    rob_g = robust_list(rs, keep=lambda t: not species_named(t))
    heat = (f'The {len(taxa)} taxa with a call at q &lt; 0.05 in any row of the table, {len(genus)} of them named only '
            "to genus or family. One column per group, laid out as the design. Each cell is the group's mean log2 share of "
            "reads over the samples where the taxon is present, which is what the abundance model fits, relative to the "
            "taxon's own average over the six groups: red is higher, blue lower. A hatched gray cell is a group with the "
            'taxon in no sample ("0/n"); "k/n" marks a group where only some samples have it. Rows are clustered by '
            "their pattern across the six groups (average linkage, Euclidean distance). The strip marks each of the 12 "
            "tests: filled for a robust call, hollow for a call lost when one mouse is left out, a dash where the taxon "
            "is missing from a group the test compares. I–C is inulin vs cellulose, R–C RS vs cellulose, R–I RS vs "
            "inulin; C, I and R are irradiation within that diet."
            + (" Every call is also a call at the stricter filter." if at6 else ""))
    forest = ("The same taxa in the same order, with every test from the table: one facet per question, its three "
              "contrasts stacked within each taxon's row. Each estimate is the log2 change in share of reads, first-named group "
              "relative to the second, with its 95% interval from MaAsLin3's t test. Filled marks are robust calls, "
              "hollow marks calls lost when one mouse is left out, faint marks no call; no mark means the taxon is "
              "missing from a group the test compares. Robust calls: " + ". ".join(rob) + ".")
    gh = (f"The same figure limited to the {len(genus)} taxa named no further than genus or family, clustered again. "
          "Nothing is re-tested.")
    gf = (f"The same tests for the {len(genus)} taxa named no further than genus or family. "
          + (("Robust calls: " + ". ".join(rob_g) + ".") if rob_g else "None has a robust call."))
    return "\n".join([
        panel(1, "Clustered factorial heatmap", f"{PRE}_factorial_heatmap",
              "Heatmap of taxa by the six Diet x Condition groups, rows clustered, with a strip marking calls in each test", heat, gh),
        panel(2, "Faceted forest plot", f"{PRE}_factorial_forest",
              "Forest plot of every test for the same taxa, one facet per question", forest, gf)])


if __name__ == "__main__":
    rs = rows()
    print(table_body(rs))
    print("<!-- panels -->")
    print(factorial_views(rs))
