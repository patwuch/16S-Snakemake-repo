#!/usr/bin/env python3
"""Check comparison-day17-day25/report.html against the tables it quotes.

Every figure the page references must exist and be non-empty, and each number
below is recomputed from its table and must appear in the report text as
formatted. Exits non-zero listing every miss.

Usage, inside the qiime2 container:
    python3 check_report_day17_day25.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-day17-day25")
TAB = os.path.join(OUT, "tables")
LN2 = np.log(2)

html = open(os.path.join(OUT, "report.html")).read()
text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
misses = []


def expect(value, why):
    if value not in text:
        misses.append(f"{why}: '{value}' not in report")


for src in re.findall(r'src="([^"]+)"', html):
    path = os.path.join(OUT, src)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        misses.append(f"figure missing or empty: {src}")

alpha = pd.read_csv(os.path.join(TAB, "diversity_day17day25_alpha_stats.tsv"), sep="\t")
alpha = alpha[alpha.db != "Silva138"]
for r in alpha.itertuples():
    expect(f"{r.n_pairs_same_direction} / 4", f"alpha {r.arm} {r.metric} mice agreeing")
    expect(f"{r.mwu_q:.2g}", f"alpha {r.arm} {r.metric} MWU q")

beta = pd.read_csv(os.path.join(TAB, "diversity_day17day25_beta_stats.tsv"), sep="\t")
beta = beta[beta.db != "Silva138"]
for r in beta.itertuples():
    expect(f"{r.R2 * 100:.0f}%", f"beta {r.arm} {r.metric} R2")
    expect(f"{r.permanova_q:.3f}", f"beta {r.arm} {r.metric} PERMANOVA q")
    expect(f"{r.permdisp_q:.2g}" if r.permdisp_q >= 0.1 else f"{r.permdisp_q:.3f}",
           f"beta {r.arm} {r.metric} PERMDISP q")
    if r.metric != "weighted_unifrac":
        expect(f"{r.mean_within_d17:.3f}", f"beta {r.arm} {r.metric} within D17")
        expect(f"{r.mean_within_d25:.3f}", f"beta {r.arm} {r.metric} within D25")


# ---- ANCOM-BC2, the McNemar tests and the turnover test are out of this report ----
for word in ("ancom", "mcnemar", "turnover"):
    if word in html.lower():
        misses.append(f"'{word}' still appears in the report")

M3 = pd.read_csv(os.path.join(TAB, "MaAsLin3_day17day25_Greengenes2.tsv"), sep="\t")
M3_ALL = pd.read_csv(os.path.join(TAB, "MaAsLin3_day17day25_Greengenes2_all_sets.tsv"), sep="\t")
NOREADS = pd.read_csv(os.path.join(TAB, "MaAsLin3_day17day25_noreads_Greengenes2.tsv"), sep="\t")
SZ = pd.read_csv(os.path.join(TAB, "DA_day17day25_Greengenes2_structural_zeros.tsv"), sep="\t")
ABUND, PREV = M3.model == "abundance", M3.model == "prevalence"


def calls(d):
    return int((d.q < 0.05).sum())


# Taxon counts in the DA opener.
n_u = M3_ALL.feature[M3_ALL.set == "U"].nunique()
n_s1 = M3_ALL.feature[M3_ALL.set == "S1"].nunique()
expect(f"which gives {n_u} taxa. {n_s1} of them are seen in each of the four arm × day groups", "U and S1 counts")
expect(f"The other {n_u - n_s1} are missing from every sample", "structural-zero count")
if SZ.taxon.nunique() != n_u - n_s1:
    misses.append("structural-zero table does not hold every taxon missing from a group")
skipped = M3[ABUND].groupby("row").skipped_median.sum()
expect(f"(up to {int(skipped.max())} taxa per comparison)", "median-test skips")
depth = pd.read_csv(os.path.join(TAB, "collapsed_counts_day17day25_Greengenes2.tsv"), sep="\t", index_col=0).sum()
expect(f"({depth.min():,}–{depth.max():,} reads per sample)", "depth range")

# MaAsLin3 question table, one html row per comparison.
LABEL = {"Day 17 to 25 in IFA": "In IFA", "Day 17 to 25 in IFA+IR": "In IFA+IR",
         "Average Day 17 to 25 shift": "Averaged over both arms", "IFA+IR vs IFA at Day 17": "At Day 17",
         "IFA+IR vs IFA at Day 25": "At Day 25", "Average IFA+IR vs IFA": "Averaged over both days",
         "Interaction (IFA+IR shift minus IFA shift)": "IFA+IR change minus IFA change"}
for row, lab in LABEL.items():
    a, p = M3[ABUND & (M3.row == row)], M3[PREV & (M3.row == row)]
    pres = "—" if p.q.notna().sum() == 0 else f"{p.q.notna().sum()} / {calls(p)}"
    qmin = M3[M3.row == row].q.min()
    qtxt = f"<b>{qmin:.3f}</b>" if qmin < 0.05 else f"{qmin:.3f}"
    cell = (f'<td>{lab}</td><td class="n">{a.q.notna().sum()} / {calls(a)}</td>'
            f'<td class="n">{pres}</td><td class="n">{qtxt}</td>')
    if cell not in html:
        misses.append(f"MaAsLin3 table row wrong: {lab}")

# Named effects.
avg = M3[ABUND & (M3.row == "Average Day 17 to 25 shift") & (M3.q < 0.05)]
expect(f"Averaged over the eight mice, {len(avg)} taxa pass", "averaged hit count")
for c in avg.coef:
    expect(f"{abs(c):.2f}", f"averaged coefficient {c:.2f}")
mur = avg[avg.feature.str.contains("f__Muribaculaceae") & (avg.coef < 0)]
expect(f"Seven Muribaculaceae fall by {mur.coef.max() * -1:.2f} to {mur.coef.min() * -1:.2f} log2", "Muribaculaceae range")
if len(mur) != 7:
    misses.append(f"text says seven Muribaculaceae fall; table has {len(mur)}")
ir = M3[ABUND & (M3.row == "Day 17 to 25 in IFA+IR")]
near = ir[ir.q == ir.q.min()]
expect(f"four taxa stop one correction step short at q = {near.q.iloc[0]:.3f}", "IFA+IR near miss")
if len(near) != 4:
    misses.append(f"text says four IFA+IR taxa at the smallest q; table has {len(near)}")
for c in near.coef:
    expect(f"{abs(c):.2f}", f"IFA+IR near-miss coefficient {c:.2f}")
ia = M3[M3.row == "Interaction (IFA+IR shift minus IFA shift)"]
expect(f"0 of {ia[ia.model == 'abundance'].q.notna().sum()} in abundance and 0 of "
       f"{ia[ia.model == 'prevalence'].q.notna().sum()} in presence, smallest q {ia.q.min():.3f}", "interaction")
if calls(ia):
    misses.append("text says the interaction finds no taxon")

# Sensitivity: without depth, and without the per-mouse term.
k = ["row", "feature", "model"]
nr = NOREADS.merge(M3, on=k, suffixes=("_nr", ""))
ir_nr = nr[(nr.row == "Day 17 to 25 in IFA+IR") & (nr.model == "abundance") & (nr.q_nr < 0.05)]
expect(f"three of the four IFA+IR taxa pass (q = {ir_nr.q_nr.max():.3f})", "no-depth IFA+IR")
if len(ir_nr) != 3 or not set(ir_nr.feature) <= set(near.feature):
    misses.append("no-depth IFA+IR calls are not three of the four near misses")
# "every coefficient above": the averaged hits and the four IFA+IR near misses.
named = nr[(nr.model == "abundance") & (((nr.row == "Average Day 17 to 25 shift") & (nr.q < 0.05))
                                        | ((nr.row == "Day 17 to 25 in IFA+IR") & nr.feature.isin(near.feature)))]
if (named.coef - named.coef_nr).abs().max() > 0.2:
    misses.append("a coefficient named in the text moves by more than 0.2 log2 without depth")
unp = {r: calls(M3[ABUND & (M3.row == r)].rename(columns={"q": "_q", "q_unpaired": "q"}))
       for r in ("Day 17 to 25 in IFA", "Day 17 to 25 in IFA+IR")}
expect(f"Without the per-mouse term, {unp['Day 17 to 25 in IFA']} taxa pass in IFA and "
       f"{unp['Day 17 to 25 in IFA+IR']} in IFA+IR", "unpaired counts")

# Largest shifts table.
cnt = pd.read_csv(os.path.join(TAB, "collapsed_counts_day17day25_Greengenes2.tsv"), sep="\t", index_col=0)
smp = pd.read_csv(os.path.join(TAB, "samples_day17day25.tsv"), sep="\t", index_col=0).loc[cnt.columns]
rel = cnt / cnt.sum() * 100
ra = {(a, d): rel.loc[:, (smp.Arm == a) & (smp.Timepoint == d)].mean(axis=1)
      for a in ("IFA", "IFA_IR") for d in ("Day17", "Day25")}
delta = sum(ra[a, "Day25"] - ra[a, "Day17"] for a in ("IFA", "IFA_IR")) / 2
delta = delta[~delta.index.str.endswith("g__;s__")]
top = delta.abs().sort_values(ascending=False).index[:11]


def pct(v):
    return "0%" if v == 0 else ("&lt;0.1%" if v < 0.05 else f"{v:.1f}%")


def qcell(row, t):
    d = M3[ABUND & (M3.row == row) & (M3.feature == t)]
    if d.empty or pd.isna(d.q.iloc[0]):
        return "—"
    q = d.q.iloc[0]
    return f"{q:.3f}" + ("†" if row.startswith("Average") and (rel.loc[t, smp.Arm == "IFA"] == 0).all() else "")


# The largest-shift table sits in Community composition, before the DA methods, so it carries
# shares only: no q-value may appear there.
_comp = html[html.index("<h2>Community composition</h2>"):html.index("</section>", html.index("<h2>Community composition</h2>"))]
_ctext = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", _comp))
for t in top:
    cells = [pct(ra["IFA", "Day17"][t]), pct(ra["IFA", "Day25"][t]), pct(ra["IFA_IR", "Day17"][t]), pct(ra["IFA_IR", "Day25"][t])]
    if " ".join(cells) not in _ctext:
        misses.append(f"composition largest-shift row {t.split(';')[-1]}: '{' '.join(cells)}' not in Community composition")
if re.search(r"\bq\b|0\.0[0-9]{2}\b", _ctext):
    misses.append("Community composition quotes a q-value; it precedes the DA methods")
# Composition bars: one toggled view per rank, Kingdom to Species, Genus shown by default,
# each view the two arms side by side.
for _rk in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"):
    _v = f"rk-{_rk[0].lower()}"
    if (f'<div class="view {_v}">' not in _comp or f'value="{_v}"' not in _comp or f".view.{_v}{{display:block}}" not in html
            or any(f"taxa_day17day25_{_a}_Greengenes2_{_rk}_groups.png" not in _comp for _a in ("IFA", "IFA_IR"))):
        misses.append(f"composition rank panel: {_rk} view missing or not wired")
if 'value="rk-g" checked' not in _comp:
    misses.append("composition rank panel: Genus is not the default view")
ms = [t for t in top if "schaedleri" in t][0]
expect(f"does not pass (averaged q {qcell('Average Day 17 to 25 shift', ms)})", "M. schaedleri averaged q")

# Structural zeros.
for _f in ("presence_day17day25_Greengenes2_heatmap.png", "presence_day17day25_Greengenes2_gplus_heatmap.png"):
    if _f not in html:
        misses.append(f"presence dot map missing: {_f}")
if "untestable_heatmap" in html:
    misses.append("the per-sample structural-zero heatmap is obsolete; use the presence dot map")
SHORT = {"Day 17 to 25 in IFA": "In IFA", "Day 17 to 25 in IFA+IR": "In IFA+IR",
         "IFA+IR vs IFA at Day 17": "At Day 17", "IFA+IR vs IFA at Day 25": "At Day 25"}
cellname = lambda c: c.replace("Day17", "Day 17").replace("Day25", "Day 25")
for con, d in SZ.groupby("contrast"):
    first, second = d.absent_in.value_counts().to_dict(), None
    row = f'<td>{SHORT[con]}</td><td class="n">{len(d)}</td>'
    if row not in html:
        misses.append(f"structural-zero count row wrong: {con}")
    for a, n in first.items():
        where = a.split()[1] if con.startswith("Day") else a.split()[0]
        where = cellname(where)
        expect(f"{n} {'at' if con.startswith('Day') else 'in'} {where}", f"structural zeros {con} missing from {a}")
    if f'<td class="n">{(d.separation == "complete").sum()}</td></tr>' not in html.split(f"<td>{SHORT[con]}</td>", 1)[1][:200]:
        misses.append(f"complete count wrong: {con}")
comp = SZ[SZ.separation == "complete"]
for r in comp.itertuples():
    expect(f"{SHORT[r.contrast]} {cellname(r.present_in)} {cellname(r.absent_in)} "
           f"{r.mean_share_where_present:.2f}% ({r.max_share:.2f}%)", f"complete row {r.label} {r.contrast}")
WORDS = {9: "nine", 13: "Thirteen"}
expect(f"{WORDS.get(len(comp), len(comp))} complete separations in {WORDS.get(comp.taxon.nunique(), comp.taxon.nunique())} taxa",
       "complete totals")
over_time = comp[comp.contrast.str.startswith("Day")]
appear = over_time.absent_in.str.endswith("Day17").sum()
expect(f"Over time, seven of the eight complete separations are taxa that appear by Day 25", "appearances")
if (len(over_time), appear) != (8, 7):
    misses.append(f"text says 7 of 8 over-time separations are appearances; table has {appear} of {len(over_time)}")
if not comp[~comp.contrast.str.startswith("Day")].absent_in.str.startswith("IFA ").all():
    misses.append("text says every between-arm complete separation is absent from IFA")
un = pd.read_csv(os.path.join(TAB, "DA_day17day25_Greengenes2_untestable_absent_from_a_cell.tsv"), sep="\t")
ab = un[[c for c in un.columns if "(n=" in c]].eq(0).rename(columns=lambda c: c.split(" (n=")[0])
only17 = ~ab["IFA Day25"] & ~ab["IFA+IR Day25"]
ifa, irr = ab["IFA Day17"], ab["IFA+IR Day17"]
expect(f"{int(only17.sum())} of the 57 structural zeros are missing only at Day 17: {int((only17 & ifa & ~irr).sum())} in the IFA "
       f"arm alone, {int((only17 & irr & ~ifa).sum())} in IFA+IR alone and {int((only17 & ifa & irr).sum())} in both",
       "heatmap caption pattern")
# Share of reads the structural-zero taxa make up together, per sample, averaged per group.
_cell = smp.Arm.map({"IFA": "IFA", "IFA_IR": "IFA+IR"}) + " " + smp.Timepoint
sh = rel.loc[un.taxon].sum().groupby(_cell).mean()
expect(f"{sh['IFA Day17']:.2f}% of reads in IFA mice at Day 17, against {sh['IFA Day25']:.2f}% at Day 25, and "
       f"{sh['IFA+IR Day17']:.2f}% and {sh['IFA+IR Day25']:.2f}% in IFA+IR", "structural-zero shares")

# Summary and DA: three taxa IFA+IR carries that are "never seen in IFA".
for _pat in ("acidifaciens", "Desulfovibrio_446353;s__Desulfovibrio sp011039135", "g__Duncaniella;s__$"):
    _t = [t for t in cnt.index if re.search(_pat, t)]
    if len(_t) != 1 or (cnt.loc[_t[0], smp.Arm == "IFA"] > 0).any() or not (cnt.loc[_t[0], smp.Arm == "IFA_IR"] > 0).any():
        misses.append(f"text says {_pat} is carried by IFA+IR and never seen in IFA; the counts disagree")
if "two views of one change" in html or "change alike</h3>" in html:
    misses.append("text closes the arm question; presence differences are open and untested")

# ---- summary section ----
if not (beta.permanova_q == beta.permanova_q.min()).all():
    misses.append("summary says every PERMANOVA row sits at the permutation floor; they differ")
expect(f"(R² {beta.R2.min() * 100:.0f}–{beta.R2.max() * 100:.0f}%)", "summary R2 range")
_spread = beta[beta.permdisp_q < 0.05].groupby("arm").metric.apply(set)
if not all(v == {"jaccard", "braycurtis", "weighted_unifrac"} for v in _spread):
    misses.append(f"summary names Jaccard, Bray-Curtis and weighted UniFrac as the spread metrics; tables give {dict(_spread)}")
expect(f"(q = {near.q.iloc[0]:.3f} against {M3[M3.row == 'Day 17 to 25 in IFA'].q.min():.3f})", "not-established q pair")

_parts = ("<h2>Differential abundance</h2>", "Structural zeros: taxa missing from a whole group",
          "A taxon missing from every sample of one group", "<th class=\"n\">Structural zeros</th>",
          "33 of the 57 structural zeros", "presence_day17day25_Greengenes2_heatmap.png", "Shares are relative abundance",
          "<th>Taxon · Greengenes2</th>", "can be compared on abundance.</p>",
          "MaAsLin3: abundance and presence", "MaAsLin3 1.5.7, on unrarefied counts", "Abundance<br>tested / hits",
          "Averaged over the eight mice", "MaAsLin3_day17day25_Greengenes2_forest_time.png",
          "Two choices cost power", 'class="finding', "<h2>Summary</h2>")
_da0 = html.find("<h2>Differential abundance</h2>")
_order = [html.find(_h, _da0) for _h in _parts]
if -1 in _order or _order != sorted(_order):
    misses.append(f"DA parts out of order or missing: {[t for t, q in zip(_parts, _order) if q == -1] or 'order'}")

# ---- beta diversity: each test's heading, explanation, table, finding, plot, then caveats ----
_b = html[html.index("<h2>Beta diversity</h2>"):html.index("</section>", html.index("<h2>Beta diversity</h2>"))]
_steps = ["PCoA: how far apart", "Points closer together", "Day 17 and Day 25 separate in every panel",
          "diversity_day17day25_beta_pcoa.png", "A PCoA is a picture, not a test", "PERMANOVA tests the separation",
          "PERMANOVA: do the timepoints differ in location?", "PERMANOVA asks whether", "<th class=\"n\">R²</th>",
          "All eight tests return p = 0.030", "diversity_day17day25_permanova_r2.png", "Do not rank the metrics",
          "PERMDISP checks\nthat", "PERMDISP: do the timepoints differ in spread?", "measures each sample's distance",
          "<th class=\"n\">within D17</th>", "Six of the eight rows show unequal spread",
          "diversity_day17day25_dispersion.png", 'class="finding']
_pos = [_b.find(t) for t in _steps]
if -1 in _pos or _pos != sorted(_pos):
    misses.append(f"beta diversity out of order or missing: {[t for t, q in zip(_steps, _pos) if q == -1] or 'order'}")

# ---- alpha diversity: explanation, table, finding, plot, caveats, findings ----
_a = html[html.index("<h2>Alpha diversity</h2>"):html.index("</section>", html.index("<h2>Alpha diversity</h2>"))]
_steps = ["Each metric is compared Day 17 against Day 25", "<th class=\"n\">MWU q</th>",
          "The three abundance-weighted metrics", "alpha_day17day25_IFA_evenness.png", "Three limits apply",
          'class="finding']
_pos = [_a.find(t) for t in _steps]
if -1 in _pos or _pos != sorted(_pos):
    misses.append(f"alpha diversity out of order or missing: {[t for t, q in zip(_steps, _pos) if q == -1] or 'order'}")

# ---- every section: its findings close it (ported from reza250909) ----
# At most four findings per section, any colour. The Summary alone may add up to two red ones
# (class "finding no": results easily misread, such as "not established") beyond its four green
# or orange; they come last, and only a heading may sit between the two groups. The design part
# is exempt: the mouse pairing needs its explanation at the very start.
SECTION_RULE_EXEMPT = {"The design"}
for _m in re.finditer(r"<section[^>]*>(.*?)</section>", html, re.S):
    _sec = _m.group(1)
    _h = re.search(r"<h2>(.*?)</h2>", _sec)
    _name = _h.group(1) if _h else "(section without a heading)"
    if _name in SECTION_RULE_EXEMPT:
        continue
    if not _sec.strip().startswith("<h2>"):
        misses.append(f"{_name}: something precedes the section heading: "
                      f"{re.sub(r'<[^>]+>', ' ', _sec.strip()[:80])!r}")
    _boxes = list(re.finditer(r'<div class="finding([^"]*)">.*?</div>', _sec, re.S))
    if not _boxes:
        continue
    _summary = _name == "Summary"
    _red = [b for b in _boxes if " no" in b.group(1)] if _summary else []
    _other = [b for b in _boxes if b not in _red]
    if len(_other) > 4:
        misses.append(f"{_name}: {len(_other)} findings, at most 4"
                      + (" green/orange" if _summary else ""))
    if len(_red) > 2:
        misses.append(f"{_name}: {len(_red)} red findings, at most 2")
    if _other and _red and _red[0].start() < _other[-1].start():
        misses.append(f"{_name}: a red finding comes before a green/orange one")
    _after = re.sub(r"\s+", " ", re.sub(r"<script.*?</script>|<[^>]+>", " ", _sec[_boxes[-1].end():], flags=re.S)).strip()
    if _after:
        misses.append(f"{_name}: content after its findings: {_after[:60]!r}")
    _rest = re.sub(r'<script.*?</script>|<div class="finding[^"]*">.*?</div>', "", _sec[_boxes[0].start():], flags=re.S)
    if _other and _red:  # a heading may introduce the red group
        _rest = re.sub(r"<h3[^>]*>.*?</h3>", "", _rest, flags=re.S)
    _between = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", _rest)).strip()
    if _between:
        misses.append(f"{_name}: findings are not together at the end ({_between[:60]!r} between them)")

# ---- Predicted function (PICRUSt2): check_picrust2.py ----
import check_picrust2  # noqa: E402
misses += check_picrust2.check(html, "day17-day25")

# LinDA is not used in this report and is not mentioned (user, 2026-10-04).
if "LinDA" in html:
    misses.append("the report mentions LinDA, which it does not use")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {len(re.findall(r'src=', html))} figures, all quoted numbers found")
