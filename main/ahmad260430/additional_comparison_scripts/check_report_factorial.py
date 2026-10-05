#!/usr/bin/env python3
"""Check comparison-Treatment-x-Radiation-factorial/report.html against its tables.

Covers the differential-abundance section and the summary lines that quote it:
every figure the page references exists and is non-empty, the ALDEx3 table's
twelve rows match aldex3_rows_factorial.py cell by cell, and each number in the
prose and captions is recomputed from its table. ANCOM-BC2 is shelved from this
report (ANALYSIS_LOG, 2026-09-30); its pre-ALDEx3 checker is archived with the
old report in ../archive_ancombc2_report_2026-09-30/. The Beta diversity
section is checked against diversity_factorial_Day17.tsv by check_beta_factorial.py
(called below); alpha diversity was audited separately (ANALYSIS_LOG, 2026-09-22)
and is not rechecked here. Exits non-zero listing every miss.

Ported from tisha251224/additional_comparison_scripts/check_report.py.

Usage, inside the qiime2 container:
    python3 check_report_factorial.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd

import aldex3_rows_factorial
from scipy.stats import fisher_exact
from plot_aldex3_factorial import genus_or_above

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Treatment-x-Radiation-factorial")
TAB = os.path.join(OUT, "tables")
STEM = "DA_factorial_Greengenes2_Day17_minpres6"

html = open(os.path.join(OUT, "report.html")).read()
text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
misses = []


def expect(value, why):
    if value not in text:
        misses.append(f"{why}: '{value}' not in report")


def cells(row_html):
    return [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", row_html)]


def fmt_row(n, up, down, q):
    u = f"+{up:.2f}" if up is not None else "—"
    d = f"−{down:.2f}" if down is not None else "—"
    return (str(n), "—" if n == 0 else f"{u} / {d}", "—" if q is None else f"{q:.1e}")


# ---- figures ----
for src in set(re.findall(r'src="([^"]+)"', html)):
    path = os.path.join(OUT, src)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        misses.append(f"figure missing or empty: {src}")
n_fig = len(set(re.findall(r'src="(figures/[^"]+\.png)"', html)))
expect(f"{n_fig} figures, PNG and SVG", "figure count")

# ---- tested set ----
main = pd.read_csv(os.path.join(TAB, f"{STEM}.tsv"), sep="\t")
strict = pd.read_csv(os.path.join(TAB, "DA_factorial_Greengenes2_Day17_minpres12.tsv"), sep="\t")
named = main.taxon.map(lambda t: bool(re.match(r"s__.", t.split(";")[-1].strip())))
expect(f"{len(main)} of them are seen in at least one mouse of every one of the six groups", "tested count")
expect(f"{named.sum()} carry a species call, and the other {(~named).sum()}", "species split")
expect(f"least 12 of the 24 mice ({len(strict)} tested)", "stricter filter taxa")
_ald = aldex3_rows_factorial.table()
if set(_ald.taxon) != set(main.taxon) or set(aldex3_rows_factorial.table(12).taxon) != set(strict.taxon):
    misses.append("ALDEx3 taxon sets differ from the tested sets")
_n_un = len(pd.read_csv(os.path.join(TAB, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t"))
expect(f"splits its power across the {len(main)} taxa seen in every group, and cannot use the {_n_un} missing from a whole group", "summary taxa count")

# ---- untestable taxa and heatmap ----
un = pd.read_csv(os.path.join(TAB, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t")
absent = un[[c for c in un.columns if "(n=" in c]].eq(0)
n_ifa = int(absent["IFA (n=4)"].sum())
expect(f"The other {len(un)} are missing from at least one whole group", "structural zero count")
# The per-mouse heatmap of the 33 was merged into the presence heatmap (peak share in each label).
if "untestable_heatmap" in html:
    misses.append("the per-mouse structural-zero heatmap is back; it is merged into the presence heatmap")
expect("Dot size is the number of mice in each group that carry it (no dot: none); the shade is the group's mean "
       "share of reads", "presence dot map encoding")
expect(f"A taxon missing from a whole group has no dot in that column; {len(un)} are, and so are left out of the abundance test",
       "presence missing taxa")
# The dot map carries no per-taxon test marks: absence shows as a missing dot.
for _w in ("◆ marks", "Grey names", "bold is the Firth"):
    if _w in html:
        misses.append(f"presence dot map caption still describes a test mark: {_w}")
expect(f"{n_ifa} of them are missing from every unirradiated IFA mouse, the emptiest column", "dot map IFA absences")
if absent.sum().idxmax() != "IFA (n=4)":
    misses.append("dot map caption: IFA is not the group missing the most taxa")

# ---- LinDA is not used in this report and is not mentioned (user, 2026-10-04); the reason ALDEx3
# was chosen is in ANALYSIS_LOG.md ----
if "LinDA" in html:
    misses.append("the report mentions LinDA, which it does not use")

# ---- ALDEx3 table ----
table = {}
question = None
i = html.index("<span style=\"text-transform:none\">γ</span> = 1</th>")
body = html[html.index("<tbody>", i):html.index("</tbody>", i)]
for r in re.findall(r"<tr[^>]*>(.*?)</tr>", body, re.S):
    c = cells(r)
    if len(c) == 8:
        question, c = c[0], c[1:]
    table[(question, c[0])] = c
A_ROWS = aldex3_rows_factorial.rows()
for r in A_ROWS:
    got = table.get((r["question"], r["label"]))
    want = fmt_row(r["n"], r["max_up"], r["max_down"], r["min_q"]) + (str(r["n_strict"]), str(r["n_gamma1"]))
    status = "none" if r["n"] == 0 else ("holds" if r["n_strict"] >= r["n"] and r["n_gamma1"] >= r["n"] else "scale-sensitive")
    if got is None:
        misses.append(f"ALDEx3 row {r['label']}: not found")
    elif tuple(got[1:6]) != want or got[6] != status:
        misses.append(f"ALDEx3 row {r['question'][:25]} {r['label']}: report {tuple(got[1:])}, tables {want + (status,)}")

# ---- ALDEx3 prose ----
_un = [r for r in A_ROWS if r["stratum"] == "unirradiated" and r["label"] == "IFA vs Control"][0]["hits"]
_ir = [r for r in A_ROWS if r["label"] == "IFA+IR vs IR"][0]["hits"]
_rad = [r for r in A_ROWS if r["row"] == "Radiation within IFA"][0]["hits"]
_int = [r for r in A_ROWS if r["row"] == "Interaction: IFA vs Control"][0]["hits"]
if not (len(_un) == len(_rad) == len(_int) == 1 and len(_ir) == 2 and
        all(h["taxon"].endswith("g__Duncaniella") for h in _un + _rad + _int) and
        [h["taxon"].split(" ")[-1] for h in _ir] == ["gordoncarteri", "freteri"]):
    misses.append("ALDEx3 hits differ from the taxa the text names")
if sum(r["n"] for r in A_ROWS) != 5 or any(r["n"] for r in A_ROWS if "Inosine" in r["label"] and "IFA" not in r["label"]):
    misses.append("text says five hits in all and none in any inosine row")
expect(f"by {_un[0]['log2FC']:.2f} log2 (q = {_un[0]['q']:.3f})", "IFA unirradiated hit")
expect(f"(+{_ir[0]['log2FC']:.2f}, q = {_ir[0]['q']:.3f})", "gordoncarteri")
expect(f"(+{_ir[1]['log2FC']:.2f}, q = {_ir[1]['q']:.3f})", "freteri")
expect(f"about {2 ** _un[0]['log2FC']:.0f}-fold", "fold-change example")
if not (9.5 <= 2 ** min(h["log2FC"] for h in _un + _ir) and 2 ** max(h["log2FC"] for h in _un + _ir) <= 15.5):
    misses.append("IFA hits are not about ten- to fifteenfold")
_g1 = aldex3_rows_factorial.table(6, "1")
_q1 = lambda row, st, t: float(aldex3_rows_factorial.pick(_g1, row, st).set_index("taxon").q[t])
_g1s = [_q1("IFA vs Control", "unirradiated", _un[0]["taxon"])] + [_q1("IFA vs Control", "irradiated", h["taxon"]) for h in _ir]
expect(f"the three treatment hits fall to q = {min(_g1s):.3f}–{max(_g1s):.3f}", "IFA hits at gamma 1")
expect(f"by {-_rad[0]['log2FC']:.2f} log2 (q = {_rad[0]['q']:.3f})", "radiation within IFA")
expect(f"(−{-_int[0]['log2FC']:.2f}, q = {_int[0]['q']:.3f})", "interaction hit")
expect(f"(q = {_q1('Radiation within IFA', '', _rad[0]['taxon']):.3f}); the interaction does not (q = "
       f"{_q1('Interaction: IFA vs Control', '', _int[0]['taxon']):.2f})", "Duncaniella at gamma 1")
_ir_rows = [r for r in A_ROWS if r["row"].startswith("Interaction")]
_nom, _tst = sum(r["nominal"] for r in _ir_rows), sum(r["tested"] for r in _ir_rows)
expect(f"{_nom} of the {_tst} interaction coefficients have raw p &lt; 0.05", "interaction nominal")
expect(f"{_nom} of {_tst} coefficients nominal where about {round(0.05 * _tst)} are expected", "summary nominal")
_mq = lambda row, st: float(aldex3_rows_factorial.pick(_ald, row, st).q.min())
expect(f"(smallest q {_mq('IFA vs Inosine', 'unirradiated'):.2f} unirradiated, {_mq('IFA vs Inosine', 'irradiated'):.2f} irradiated)",
       "IFA vs Inosine smallest q")
_iu = aldex3_rows_factorial.pick(_ald, "Inosine vs Control", "unirradiated").sort_values("p").iloc[0]
_ii = aldex3_rows_factorial.pick(_ald, "Inosine vs Control", "irradiated").sort_values("p").iloc[0]
if not (_iu.taxon.endswith("g__Alistipes_A_871400;s__") and _ii.taxon.endswith("Duncaniella freteri")):
    misses.append("inosine lead taxa differ from the text")
expect(f"{-_iu.log2FC:.2f} log2 lower than Control in unirradiated mice (raw p = {_iu.p:.3f}, q = {_iu.q:.2f})", "inosine unirradiated")
expect(f"{_ii.log2FC:.2f} log2 higher under inosine (q = {_ii.q:.3f})", "inosine irradiated")

# ---- claims added in the 2026-10-01 wording pass ----
_gm = pd.read_csv(os.path.join(TAB, "ALDEx3_factorial_Greengenes2_Day17_minpres6_gamma0.5_groupmeans.tsv"), sep="\t").set_index("taxon")
_dg = _gm.loc[_un[0]["taxon"]]
if not (_dg["IFA"] == _dg.max() and _dg["IFA+IR"] < _dg["Control"]):
    misses.append("text says Duncaniella (genus only) peaks in IFA and ends below Control in IFA+IR")
_mur = _gm.loc[[t for t in _gm.index if "f__Muribaculaceae" in t]]
_mc = _mur.sub(_mur.mean(axis=1), axis=0).mean()
if not (_mc["IFA"] == _mc.max() and 0 < _mc["IFA+IR"] < _mc["IFA"]):
    misses.append("text says the Muribaculaceae are highest in unirradiated IFA and less raised in IFA+IR")
expect(f"the taxon is missing from {int(_iu.n_zero_mice)} of the 24 mice", "Alistipes zero count")
_hits = lambda d: {(r.row, r.stratum, r.taxon) for r in d[d.q < 0.05].itertuples()}
_h5 = _hits(_ald)
_h25, _h0 = _hits(aldex3_rows_factorial.table(6, "0.25")), _hits(aldex3_rows_factorial.table(6, "0"))
_extra = _h0 - _h5
if not (_h25 == _h5 and not (_h5 - _h0) and len(_extra) == 1 and
        next(iter(_extra))[:2] == ("Inosine vs Control", "irradiated") and next(iter(_extra))[2].endswith("Duncaniella freteri")):
    misses.append("text says gamma 0.25 gives the same hits and gamma 0 adds D. freteri in Inosine+IR vs IR")
if re.search(r"log of zero|genus-level, no species call", html):
    misses.append("stale wording: 'log of zero' or 'genus-level, no species call'")

# ---- reading guide ("Both together") ----
_gi = html.find("Reading the two kinds of result together")
if not (html.find("Both together") < _gi < html.find("Which taxa move with which factor?")):
    misses.append("reading guide missing or not at the start of Both together")
_guide = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html[_gi:html.find("Which taxa move with which factor?")]))
for _w, _why in ((f"All {len(main) + len(un)}, in the dot map; the {len(un)} missing from a whole group are the empty cells", "guide taxa"),
                 (f"The {len(main)} seen in every group", "guide tested"),
                 (f"a complete separation reaches p = {fisher_exact([[4, 0], [0, 4]])[1]:.3f} whatever its size", "guide floor"),
                 (f"the strongest hit here has q = {min(h['q'] for r in A_ROWS for h in r['hits']):.3f}", "guide min q"),
                 (f"{n_ifa} of the {len(un)} taxa missing from a whole group are absent from every unirradiated IFA mouse", "guide IFA absent"),
                 (f"ALDEx3's {sum(r['n'] for r in A_ROWS)} hits", "guide hit count"),
                 ("while absences are shown, not tested", "guide absences")):
    if _w not in _guide:
        misses.append(f"{_why}: '{_w}' not in the reading guide")

# ---- figures and their captions ----
_n_g = int(sum(genus_or_above(t) for t in main.taxon))
expect(f"limited to the {_n_g} taxa named no further than genus", "genus forest count")
expect(f"The same heatmap for the {_n_g} taxa named", "genus heatmap count")
_cnt = pd.read_csv(os.path.join(TAB, "collapsed_counts_Greengenes2_Day17.tsv"), sep="\t", index_col=0)
_79 = _cnt[(_cnt > 0).sum(axis=1) >= 6]
if len(_79) != len(main) + len(un):
    misses.append(f"presence heatmap has {len(_79)} taxa, tested + untestable is {len(main) + len(un)}")
expect(f"The same, for the {int(sum(genus_or_above(t) for t in _79.index))} taxa named", "genus presence count")
for _k in ("treatment_unirradiated", "treatment_irradiated", "radiation", "interaction"):
    for _v in ("", "_gplus"):
        if f"figures/ALDEx3_factorial_Greengenes2_Day17_minpres6_gamma0.5{_v}_forest_{_k}.png" not in html:
            misses.append(f"forest plot {_k}{_v} not shown")
for _f in ("ALDEx3_factorial_Greengenes2_Day17_minpres6_gamma0.5_heatmap", "ALDEx3_factorial_Greengenes2_Day17_minpres6_gamma0.5_gplus_heatmap",
           "presence_factorial_Greengenes2_Day17_minpres6_heatmap", "presence_factorial_Greengenes2_Day17_minpres6_gplus_heatmap"):
    if f"figures/{_f}.png" not in html:
        misses.append(f"{_f} not shown")
# Data before model: the presence dot map leads, the ALDEx3 heatmap follows.
if not (0 <= html.find("figures/presence_factorial_Greengenes2_Day17_minpres6_heatmap.png")
        < html.find("figures/ALDEx3_factorial_Greengenes2_Day17_minpres6_gamma0.5_heatmap.png")):
    misses.append("the presence dot map should come before the ALDEx3 heatmap")
# ANCOM-BC2 is out of this report entirely: no result, figure or file path may name it.
for _m in re.finditer(r".{0,60}ancom.{0,60}", html, re.I | re.S):
    misses.append(f"ANCOM-BC2 still mentioned: {' '.join(_m.group(0).split())}")

# ---- structural zeros: the dot map only (2026-10-01) ----
# Taxa missing from a whole group are shown, not tested; the section is the dot map and nothing else.
expect(f"reaches the same smallest p-value ({fisher_exact([[4, 0], [0, 4]])[1]:.3f}) whatever its abundance", "design floor")
# Scoped to DA onward: Beta diversity's PERMANOVA legitimately speaks of shuffled labels and permutation p.
for _m in re.finditer(r".{0,50}(?:Fisher|Firth|shuffled labels|permutation p|2\.6 expected|q = 0\.0037).{0,50}",
                      html[html.index("<h2>Differential abundance</h2>"):], re.S):
    misses.append(f"per-taxon structural-zero test still in the report: {' '.join(_m.group(0).split())}")
# The DA-subset binary-Jaccard PERMANOVA is dropped; the interaction's community evidence is Beta diversity.
for _m in re.finditer(r".{0,50}(?:binary[- ]Jaccard|0\.0057|R² 13%).{0,50}", html, re.S):
    misses.append(f"dropped presence PERMANOVA still cited: {' '.join(_m.group(0).split())}")
_szs = html[html.index('<h3>Structural zeros: taxa present'):html.index('<h3>ALDEx3: abundance in the 46 taxa')]
# Subsections are plain h3; a figure's question sits inside its panel, never as a heading of its own.
_dahtml = html[html.index("<h2>Differential abundance</h2>"):html.index("<h2>Summary</h2>")]
if 'class="figq"' in _dahtml:
    misses.append("the DA section has a figure-question heading styled like a subsection")
for _q in ("Which taxa move with which factor?", "Do the treatments differ in unirradiated mice?",
           "Do the treatments differ in irradiated mice?", "Does radiation change the microbiome within each treatment?",
           "Does radiation's effect differ between treatments?"):
    if f'<div class="panel-q">{_q}</div>' not in _dahtml:
        misses.append(f"figure question not inside its panel: {_q}")
_szs_out = re.sub(r"<figure.*?</figure>", "", _szs, flags=re.S)
if _szs.count("<figure") != 1 or "presence_factorial_Greengenes2_Day17_minpres6_heatmap.png" not in _szs \
        or re.search(r"<p>|<table", _szs_out):
    misses.append("the structural-zero section should hold the presence dot map and nothing else")
# Presence patterns the findings and summary read off the dot map.
_dep = pd.read_csv(os.path.join(TAB, f"{STEM}_structural_zeros_depth.tsv"), sep="\t")
_un = pd.read_csv(os.path.join(TAB, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t").set_index("taxon")
_cc = [c for c in _un.columns if "(n=" in c]
_c = lambda k: _un.loc[[t for t in _un.index if t.endswith(k)][0], _cc].rename(lambda c: c.split(" (")[0])
for _k in ("Alloprevotella sp002933955", "Desulfovibrio sp011039135", "g__Duncaniella;s__"):
    if not (_c(_k)["IR"] == 4 and _c(_k)["Control"] == 0):
        misses.append(f"{_k}: not in every IR mouse and no Control mouse")
_dn = _c("g__Duncaniella;s__")
if tuple(_dn[["Control", "IR", "IFA", "IFA+IR", "Inosine", "Inosine+IR"]]) != (0, 4, 0, 4, 4, 0):
    misses.append(f"Duncaniella presence is not 0->4 in Control and IFA, 4->0 in Inosine: {dict(_dn)}")
_ba = _c("acidifaciens")
if tuple(_ba[["Control", "IR", "IFA", "IFA+IR", "Inosine", "Inosine+IR"]]) != (4, 0, 0, 3, 0, 0):
    misses.append(f"B. acidifaciens presence differs from the text: {dict(_ba)}")
_sh = lambda t: t.split(";")[-1] if not t.endswith("__") else t.split(";")[-2]
_ino = {(_sh(r.taxon), r.present_in.replace("Inosine", "X")) for r in _dep[_dep.contrast == "Inosine vs Control"].itertuples()}
_irc = {(_sh(r.taxon), r.present_in.replace("IR", "X")) for r in _dep[_dep.contrast == "Control"].itertuples()}
if len(_ino) != 4 or not _ino <= _irc:
    misses.append("text says inosine switches four of irradiation's Control taxa, same direction")
if {t for t, _ in _ino} != {"s__Bacteroides_H_857956 acidifaciens", "s__Alloprevotella sp002933955", "g__Duncaniella",
                            "s__Bacteroides_H_857956 caccae"}:
    misses.append(f"finding 3 names the wrong four inosine taxa: {sorted(t for t, _ in _ino)}")
_in_every = {t for t, where in _ino if where == "X"}
if len(_in_every) != 3:
    misses.append("finding 3: three of the four should be in every Inosine mouse, one in none")
expect("Inosine without irradiation does the same to four of them", "finding 3 inosine")
expect("The interaction, asked six ways", "six ways")
for _h in ("How taxa are tested", "Structural zeros: taxa present in some groups", "ALDEx3: abundance in the 46 taxa", "Both together"):
    if _h not in html:
        misses.append(f"DA part heading missing: {_h}")
_da0 = html.index("<h2>Differential abundance</h2>")
_order = [html.find(_h, _da0) for _h in ("How taxa are tested", "Structural zeros: taxa present", "ALDEx3: abundance in", ">Both together</h3>", "<h3>Findings</h3>", "<h2>Summary</h2>")]
if -1 in _order:
    misses.append("a DA part heading is missing")
_dsec = html[_da0:html.index("<h2>Summary</h2>")]
if 'class="note"' in _dsec:
    misses.append("the DA section has a callout box; results there are plain text")
# ---- community composition: group means of the species that differ most, shares only ----
import biom  # noqa: E402
_meta = pd.read_csv(os.path.join(TAB, "metadata_day17.tsv"), sep="\t", index_col=0)
_cnt = biom.load_table(os.path.join(os.path.dirname(OUT), "tables", "study-seqs.biom")).to_dataframe(dense=True)
_tax = pd.read_csv(os.path.join(os.path.dirname(OUT), "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                   sep="\t", index_col=0)["Taxon"]
_c = _cnt.loc[_cnt.index.intersection(_tax.index)]
_sp = _c.groupby(_tax.loc[_c.index].map(lambda t: ";".join([p for p in (x.strip() for x in t.split(";")) if p][:7]))).sum()[_meta.index]
_rel = _sp / _sp.sum() * 100
_G = ["Control", "IFA", "Inosine", "IR", "IFA+IR", "Inosine+IR"]
_mean = pd.DataFrame({g: _rel.loc[:, _meta.Group == g].mean(axis=1) for g in _G})
_mean = _mean[~_mean.index.str.endswith("g__;s__")]
_top = (_mean.max(axis=1) - _mean.min(axis=1)).sort_values(ascending=False).index[:12]
_c0 = html.index("<h2>Community composition</h2>")
_cs = html[_c0:html.index("</section>", _c0)]
_rows = re.findall(r"<tr><td>.*?</tr>", _cs)
_pct = lambda v: "&lt;0.1%" if v < 0.05 else f"{v:.1f}%"
if len(_rows) != 12:
    misses.append(f"composition table has {len(_rows)} rows, expected 12")
for _r, _t in zip(_rows, _top):
    _got = re.findall(r'<td class="n">(.*?)</td>', _r)
    if _got != [_pct(v) for v in _mean.loc[_t, _G]]:
        misses.append(f"composition row {_t.split(';')[-1]}: report {_got}, tables {[_pct(v) for v in _mean.loc[_t, _G]]}")
_s = lambda k: _mean.loc[[t for t in _mean.index if t.endswith(k)][0]]
_rng = lambda v: f"{v.min():.1f}–{v.max():.1f}%"
_dm, _mg, _dd, _ba = _s("Duncaniella muricolitica"), _s("Muribaculum gordoncarteri"), _s("Duncaniella dubosii"), _s("acidifaciens")
_ifa = ["IFA", "IFA+IR"]
expect(f"{_dm['IFA']:.1f}% of reads in IFA mice and {_dm['IFA+IR']:.1f}% in irradiated IFA mice, against {_rng(_dm.drop(_ifa))} in every other group", "composition D. muricolitica")
expect(f"{_dm['IFA']:.1f}% and {_dm['IFA+IR']:.1f}% of reads in the two IFA groups, against {_rng(_dm.drop(_ifa))} elsewhere", "summary composition")
expect(f"{_mg['IFA']:.1f}% and {_mg['IFA+IR']:.1f}% against {_rng(_mg.drop(_ifa))}", "composition M. gordoncarteri")
expect(f"({_dd['IFA']:.1f}%, against {_rng(_dd.drop('IFA'))})", "composition D. dubosii")
expect(f"only in control mice ({_ba['Control']:.1f}%)", "composition B. acidifaciens")
if _ba.drop("Control").max() >= 1:
    misses.append("composition: B. acidifaciens is not control-only")
_ph = _c.groupby(_tax.loc[_c.index].map(lambda t: t.split(";")[1].strip() if ";" in t else "unassigned")).sum()[_meta.index]
_ph = _ph / _ph.sum() * 100
_bd = pd.Series({g: _ph.loc["p__Bacteroidota", _meta.index[_meta.Group == g]].mean() for g in _G})
expect(f"makes up {_bd['IFA']:.1f}% of reads in IFA mice, against {_bd['Control']:.1f}% in control mice and {_bd['IR']:.1f}% in irradiated mice", "composition Bacteroidota")
if _bd.idxmax() != "IFA" or _bd.idxmin() != "IR":
    misses.append("composition: Bacteroidota is not highest in IFA and lowest in IR")
if re.search(r"\bq\b|\bq =|raw p", re.sub(r"<[^>]+>", " ", _cs)):
    misses.append("Community composition quotes a test result; it precedes the DA methods")
for _rk in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"):
    _v = f"rk-{_rk[0].lower()}"
    if (f'<div class="view {_v}">' not in _cs or f"taxa_factorial_Greengenes2_{_rk}_groups.png" not in _cs
            or f'value="{_v}"' not in _cs or f".view.{_v}{{display:block}}" not in html):
        misses.append(f"composition rank panel: {_rk} view missing or not wired")
if 'value="rk-g" checked' not in _cs:
    misses.append("composition rank panel: Genus is not the default view")
_ord = [html.find(x) for x in ("<h2>Beta diversity</h2>", "<h2>Community composition</h2>", "<th>Species · Greengenes2</th>",
                               "The twelve species whose", '<figure class="panel ranks">', "per-mouse bars", "<h2>Differential abundance</h2>")]
if -1 in _ord or _ord != sorted(_ord):
    misses.append(f"composition parts missing or out of order: {_ord}")

# Findings only ever at the very end of a section: after a section's first finding, nothing but
# findings (and a "Findings" heading before them). The DA section closes on at most four.
for _sec in re.split(r"(?=<section)", html):
    _t = re.search(r"<h2>(.*?)</h2>", _sec)
    if not _t or '<div class="finding' not in _sec:
        continue
    _rest = _sec[_sec.index('<div class="finding'):]
    _rest = re.sub(r'<div class="finding[^"]*">.*?</div>', "", _rest, flags=re.S)
    if re.search(r"<figure|<table|<h3|<p[ >]|<div class=\"col", _rest):
        misses.append(f"{_t.group(1)}: content follows a finding; findings must close the section")
    _nf = _sec.count('<div class="finding')
    if _t.group(1) == "Differential abundance" and _nf > 4:
        misses.append(f"Differential abundance has {_nf} findings; at most 4")
if _order != sorted(_order):
    misses.append("DA parts are out of order")

# ---- Beta diversity: tables cell by cell, and every number the beta prose quotes ----
import check_beta_factorial  # noqa: E402

misses += check_beta_factorial.check(html)

# ---- Predicted function (PICRUSt2): check_picrust2.py ----
import check_picrust2  # noqa: E402
misses += check_picrust2.check(html, "factorial")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {n_fig} figures, all quoted numbers found")
