#!/usr/bin/env python3
"""The PICRUSt2 section of both ahmad260430 reports, built from the tables.

Everything the section says is computed here from the PICRUSt2 output (picrust2/)
and the tables written by run_picrust2_linda.R, picrust2_permanova.R,
plot_picrust2.R and picrust2_acetate.py. Each qualitative claim the prose makes
("below chance", "holds without mouse X") is asserted against the numbers
before the section is rendered, so a rerun that changes a claim fails loudly
instead of leaving stale prose behind.

    python3 picrust2_report.py factorial|day17-day25 --write
        (re)writes the section and its summary, interaction and footer lines
        into the report, between marker comments; running it twice is a no-op
    python3 picrust2_report.py factorial|day17-day25
        exits non-zero if the report's PICRUSt2 parts differ from what the
        tables give (the check each report's checker runs)

Inside the qiime2 container.
"""
import html as _html
import os
import re
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
PIC = os.path.join(STUDY_DIR, "picrust2", "pipeline")
LEVELS = ("MetaCyc", "KO", "EC")
VIEW = {"MetaCyc": "mc", "KO": "ko", "EC": "ec"}
LEVEL_NAME = {"MetaCyc": "MetaCyc pathways", "KO": "KEGG orthologs", "EC": "enzymes"}
LEVEL_NOUN = {"MetaCyc": "pathways", "KO": "orthologs", "EC": "enzymes"}
FILES = {"MetaCyc": "pathways_out/path_abun_unstrat.tsv.gz",
         "KO": "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
         "EC": "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"}


@dataclass(frozen=True)
class Contrast:
    name: str       # as in the LinDA tables, first-named group positive
    question: str   # the report's question this contrast answers
    label: str      # its row label under that question


@dataclass(frozen=True)
class Report:
    key: str
    dir: str
    contrasts: tuple
    sweep: tuple    # prevalence filters refitted, 0 = every feature
    top: int        # rows per contrast in the results tables
    heatmap_top: int

    @property
    def out(self):
        return os.path.join(STUDY_DIR, self.dir)

    @property
    def half(self):
        return max(self.sweep)


FACTORIAL = Report(
    key="factorial", dir="comparison-Treatment-x-Radiation-factorial",
    contrasts=(
        Contrast("IFA vs Control", "Do the treatments differ in unirradiated mice?", "IFA vs Control"),
        Contrast("Inosine vs Control", "Do the treatments differ in unirradiated mice?", "Inosine vs Control"),
        Contrast("IFA vs Inosine", "Do the treatments differ in unirradiated mice?", "IFA vs Inosine"),
        Contrast("IFA+IR vs IR", "Do the treatments differ in irradiated mice?", "IFA+IR vs IR"),
        Contrast("Inosine+IR vs IR", "Do the treatments differ in irradiated mice?", "Inosine+IR vs IR"),
        Contrast("IFA+IR vs Inosine+IR", "Do the treatments differ in irradiated mice?", "IFA+IR vs Inosine+IR"),
        Contrast("IR vs Control", "Does radiation change the microbiome within each treatment?", "Control"),
        Contrast("IFA+IR vs IFA", "Does radiation change the microbiome within each treatment?", "IFA"),
        Contrast("Inosine+IR vs Inosine", "Does radiation change the microbiome within each treatment?", "Inosine"),
    ),
    sweep=(0, 6, 9, 12), top=3, heatmap_top=5)
DAY = Report(
    key="day17-day25", dir="comparison-day17-day25",
    contrasts=(
        Contrast("IFA Day25 vs IFA Day17", "Does the community change from Day 17 to Day 25 within each arm?", "IFA"),
        Contrast("IFA+IR Day25 vs IFA+IR Day17", "Does the community change from Day 17 to Day 25 within each arm?", "IFA+IR"),
        Contrast("Shift, IFA+IR vs IFA", "Does the shift differ between arms?", "IFA+IR vs IFA"),
    ),
    sweep=(0, 4, 6, 8), top=5, heatmap_top=10)
REPORTS = {r.key: r for r in (FACTORIAL, DAY)}


def esc(s):
    return _html.escape(str(s), quote=False)


def fmtp(x):
    """p or q as the reports print it; a fourth decimal where three would round across 0.05."""
    if x < 0.001:
        return f"{x:.1e}"
    if x < 0.1:
        return f"{x:.4f}" if (x < 0.05) != (round(x, 3) < 0.05) else f"{x:.3f}"
    return f"{x:.2f}"


def chance(n):
    return round(0.05 * n)


def rng(vals, fmt):
    lo, hi = min(vals), max(vals)
    return f"{lo:{fmt}}" if f"{lo:{fmt}}" == f"{hi:{fmt}}" else f"{lo:{fmt}}–{hi:{fmt}}"


def read_linda(r, lv, suffix=""):
    return pd.read_csv(os.path.join(r.out, "tables", f"PICRUSt2_{lv}_linda{suffix}.tsv"), sep="\t", quoting=3)


def minpres(mp):
    return f"_minpres{mp}" if mp else ""


def load(r):
    """Every fit the section quotes, keyed (level, min presence, exclusion suffix)."""
    fits = {}
    for lv in LEVELS:
        for mp in r.sweep:
            fits[(lv, mp, "")] = read_linda(r, lv, minpres(mp))
        for f in os.listdir(os.path.join(r.out, "tables")):
            m = re.fullmatch(rf"PICRUSt2_{lv}_linda_minpres{r.half}(_excl.+)\.tsv", f)
            if m:
                fits[(lv, r.half, m.group(1))] = read_linda(r, lv, f"_minpres{r.half}{m.group(1)}")
    return fits


def validate(r, fits):
    """Properties every fit must have before anything is quoted from it."""
    bad = []
    for (lv, mp, ex), d in fits.items():
        if set(d.contrast) != {c.name for c in r.contrasts}:
            bad.append(f"{lv} {mp}{ex}: contrasts {sorted(set(d.contrast))}")
        for c in r.contrasts:
            x = d[d.contrast == c.name]
            if len(x) and not np.allclose(multipletests(x.p, method="fdr_bh")[1], x.q):
                bad.append(f"{lv} {mp}{ex} {c.name}: q is not BH within the contrast")
        if mp and (d.n_predicted < mp).any():
            bad.append(f"{lv} {mp}{ex}: holds features predicted in fewer than {mp} samples")
        if d.description.isna().any():
            bad.append(f"{lv} {mp}{ex}: features without a description")
    return bad


def stats(d, contrast):
    x = d[d.contrast == contrast]
    s = x[x.p < 0.05]
    return dict(n=len(x), raw=len(s), down=int((s.log2fc < 0).sum()), up=int((s.log2fc > 0).sum()),
                q=int((x.q < 0.05).sum()), med=float(x.log2fc.median()))



def samples(r):
    """Sample -> group, as the LinDA fits group them."""
    if r.key == "factorial":
        m = pd.read_csv(os.path.join(STUDY_DIR, "metadata.tsv"), sep="\t", index_col=0, keep_default_na=False)
        return m[m.Timepoint == "Day17"].Group
    return pd.read_csv(os.path.join(r.out, "metadata_joint_cells.tsv"), sep="\t", index_col=0, keep_default_na=False).Group


def sides(contrast):
    """The two groups a contrast compares; the arm shift compares the arms at Day 25."""
    g = contrast.removeprefix("Shift, ").split(" vs ")
    return [x if not contrast.startswith("Shift, ") else f"{x} Day25" for x in g]


_REL = {}


def rel(lv):
    if lv not in _REL:
        ab = pd.read_csv(os.path.join(PIC, FILES[lv]), sep="\t", index_col=0)
        _REL[lv] = ab / ab.sum()
    return _REL[lv]


def dominant(r, fits, contrast, lv="KO"):
    """The mouse highest on most of a contrast's flagged features (raw p < 0.05, fit among
    features predicted in at least half the samples), on the side with more of them.
    Returns (mouse, features it tops, flagged features on that side, side)."""
    grp = samples(r)
    d = fits[(lv, r.half, "")]
    d = d[(d.contrast == contrast) & (d.p < 0.05)]
    g1, g2 = sides(contrast)
    up, down = d[d.log2fc > 0].feature, d[d.log2fc < 0].feature
    side, feats = (g1, up) if len(up) >= len(down) else (g2, down)
    top = rel(lv).loc[feats, grp.index[grp == side]].idxmax(axis=1).value_counts()
    return top.index[0], int(top.iloc[0]), len(feats), side

def seg(name, first_title):
    radios = "".join(f'<input type="radio" name="{name}" id="{name}{VIEW[lv]}" value="{VIEW[lv]}"{" checked" if lv == "MetaCyc" else ""}>'
                     f'<label for="{name}{VIEW[lv]}">{lv}</label>' for lv in LEVELS)
    return (f'<figure class="panel">\n  <div class="panel-head"><span class="panel-title">{first_title}</span>\n'
            f'    <span class="seg">{radios}</span></div>\n')


def table(head, rows, cls="long pic"):
    th = "".join(('<th class="n">' if n else "<th>") + f"{h}</th>" for h, n in head)
    body = "\n".join("<tr>" + "".join(c if c.startswith("<td") else f"<td>{c}</td>" for c in row) + "</tr>" for row in rows)
    return f'<div class="scroll">\n<table class="{cls}">\n<thead><tr>{th}</tr></thead>\n<tbody>\n{body}\n</tbody></table>\n</div>\n'


def num(x):
    return f'<td class="n">{x}</td>'


def status(qs):
    if all(qs):
        return '<span class="chip ok">q &lt; 0.05 at every filter</span>'
    return '<span class="chip warn">filter-dependent</span>' if any(qs) else '<span class="chip warn">below threshold</span>'


def feature_table(r, fits):
    head = [("Database", False), ("Contrast", False), ("raw p &lt; 0.05", True),
            (f"raw p &lt; 0.05, ≥ {r.half} samples (chance)", True),
            ("q &lt; 0.05 by filter: " + " / ".join(str(m) if m else "all" for m in r.sweep), True), ("Status", False)]
    rows = []
    for lv in LEVELS:
        for i, c in enumerate(r.contrasts):
            a, w = stats(fits[(lv, 0, "")], c.name), stats(fits[(lv, r.half, "")], c.name)
            qs = [stats(fits[(lv, m, "")], c.name)["q"] for m in r.sweep]
            first = [f'<td rowspan="{len(r.contrasts)}">{LEVEL_NAME[lv][0].upper() + LEVEL_NAME[lv][1:]}<br><span class="sub">{a["n"]:,} features</span></td>'] if i == 0 else []
            rows.append(first + [esc(c.name), num(f"{a['raw']:,}"), num(f"{w['raw']:,} / {w['n']:,} (~{chance(w['n']):,})"),
                                 num(" / ".join(f"{q:,}" for q in qs)), status(qs)])
    return table(head, rows)


def description(lv, d):
    return "no name in either KO map" if d == "not_found" else d


def top_table(r, fits, lv):
    head = [("Contrast", False), ("ID", False), ("Description", False), ("log2FC", True), ("raw p", True),
            (f"q, ≥ {r.half} samples", True), ("q, all features", True)]
    allf = fits[(lv, 0, "")].set_index(["contrast", "feature"])
    rows = []
    for c in r.contrasts:
        x = fits[(lv, r.half, "")]
        x = x[x.contrast == c.name].sort_values("p").head(r.top)
        for i, t in enumerate(x.itertuples()):
            qa = fmtp(allf.loc[(c.name, t.feature)].q) if (c.name, t.feature) in allf.index else "—"
            first = [f'<td rowspan="{len(x)}">{esc(c.name)}</td>'] if i == 0 else []
            rows.append(first + [f"<code>{esc(t.feature.replace('ko:', ''))}</code>", esc(description(lv, t.description)),
                                 num(f"{t.log2fc:+.2f}"), num(fmtp(t.p)), num(fmtp(t.q)), num(qa)])
    return table(head, rows)


def figure_panel(name, title, stem, caption):
    out = seg(name, title)
    for lv in LEVELS:
        out += (f'    <div class="view {VIEW[lv]}">\n      <div class="one">\n'
                f'      <img src="figures/picrust2_{lv}_{stem}.png" alt="{esc(title)} ({lv})">\n      </div>\n'
                f'      <figcaption>{caption(lv)}</figcaption>\n    </div>\n')
    return out + "</figure>\n"


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in REPORTS:
        sys.exit(f"usage: picrust2_report.py {'|'.join(REPORTS)} [--write]")
    rep = REPORTS[sys.argv[1]]
    fits = load(rep)
    bad = validate(rep, fits)
    if bad:
        sys.exit("\n".join(bad))
    if "--dominant" in sys.argv:
        for c in rep.contrasts:
            print(c.name, *dominant(rep, fits, c.name), sep="\t")
