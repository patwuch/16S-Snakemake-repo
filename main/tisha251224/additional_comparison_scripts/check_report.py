#!/usr/bin/env python3
"""Check comparison-Diet-x-Condition-factorial/report.html against its tables.

Every figure the page references must exist and be non-empty, and each number
below is recomputed from its table and must appear in the report text as
formatted. Exits non-zero listing every miss.

Usage, inside the qiime2 container:
    python3 check_report.py [report.html]
"""
import os
import re
import sys

import numpy as np
import pandas as pd

from da_rows import rows, species_named
from plot_interaction_forest import PANELS, joint_q, load_tables

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
TAB = os.path.join(OUT, "tables")

# An optional path checks another copy (mutation tests); the default is the report.
html = open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(OUT, "report.html")).read()
text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
misses = []


def expect(value, why):
    if value not in text:
        misses.append(f"{why}: '{value}' not in report")


for src in set(re.findall(r'src="([^"]+)"', html)):
    path = os.path.join(OUT, src)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        misses.append(f"figure missing or empty: {src}")
n_fig = len(set(re.findall(r'src="(figures/[^"]+\.png)"', html)))
expect(f"{n_fig} figures, PNG and SVG", "figure count")

# ---- design ----
meta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t")
expect(f"Samples {len(meta)}", "sample count")

# ---- diversity ----
div = pd.read_csv(os.path.join(TAB, "diversity_factorial.tsv"), sep="\t")
alpha = div[div.analysis == "alpha"]
for r in alpha.itertuples():
    q = r.p_BH
    expect(f"{q:.1e}" if q < 1e-3 else f"{q:.2g}", f"alpha {r.metric} {r.term} q")

perm = div[div.analysis == "beta_permanova"].set_index(["metric", "term"])
for (metric, term), r in perm.iterrows():
    if term == "Residual":
        continue
    q = r.p_BH
    expect("0.0001" if q <= 1e-4 else (f"{q:.4f}" if abs(q - 0.05) < 0.001 else f"{q:.2g}"),
           f"PERMANOVA {metric} {term} q")
    expect(f"{r.effect * 100:.0f}%", f"PERMANOVA {metric} {term} R2")
lo = lambda t: perm.xs(t, level="term").effect
expect(f"{lo('Diet').min() * 100:.0f}–{lo('Diet').max() * 100:.0f}%", "Diet R2 range")
expect(f"{lo('Residual').min() * 100:.0f}–{lo('Residual').max() * 100:.0f}%", "Residual R2 range")

for a in ("beta_dispersion_factorial",):
    for r in div[div.analysis == a].itertuples():
        q = r.p_BH
        expect(f"{q:.4f}" if q < 0.01 else (f"{q:.3f}" if q < 0.1 else f"{q:.2f}"), f"{a} {r.metric} {r.term} q")
for r in div[(div.analysis == "beta_dispersion") & (div.term == "Cell")].itertuples():
    q = r.p_BH
    expect(f"{q:.4f}" if q < 0.01 else f"{q:.2f}", f"PERMDISP cell {r.metric} q")

# ---- differential abundance: taxon sets ----
PRES = "presence_factorial_Greengenes2"
tested_every = pd.read_csv(os.path.join(TAB, f"{PRES}.tsv"), sep="\t")
untestable = pd.read_csv(os.path.join(TAB, f"{PRES}_untestable_absent_from_a_cell.tsv"), sep="\t")
ma_counts = pd.read_csv(os.path.join(TAB, "maaslin3", "counts_Greengenes2.tsv"), sep="\t", index_col=0)
n2, n6 = int(((ma_counts > 0).sum() >= 2).sum()), int(((ma_counts > 0).sum() >= 6).sum())
if n6 != len(tested_every) + len(untestable):
    misses.append(f"presence set: {n6} taxa in >= 6 samples, lists give {len(tested_every) + len(untestable)}")
expect(f"Of 188 features, {n2} are seen in at least two of the 19 samples and are tested for abundance with MaAsLin3", "MaAsLin3 taxon set")
expect(f"{n6} are seen in at least six of the 19 samples", "presence taxon set")
expect(f"The other {len(untestable)} of those {n6} are missing from at least one whole group", "structural-zero count")
expect(f"Of the {len(untestable)} structural zeros", "untestable count")
rsir = untestable.columns[untestable.columns.str.startswith("RS+IR")][0]
rs = untestable.columns[untestable.columns.str.startswith("RS (")][0]
absent = untestable[untestable[rsir] == 0]
expect(f"{len(absent)} of the {len(untestable)} are absent from all three RS+IR", "absent from RS+IR")
expect(f"{int((absent[rs] > 0).sum())} of those {len(absent)} are present in unirradiated RS", "present in RS")

# ---- MaAsLin3 table: every cell and chip from run_maaslin3_rows.R's tables ----
MINUS = "\u2212"
num = lambda v: f"+{v:.2f}" if v >= 0 else f"{MINUS}{-v:.2f}"
_da0 = html.index("<th class=\"n\">Typical taxon</th>")
_da_html = html[_da0:html.index("</table>", _da0)]
da = rows()
_trs = re.findall(r"<tr[^>]*>(.*?)</tr>", _da_html, re.S)
if len(_trs) != len(da):
    misses.append(f"MaAsLin3 table: {len(_trs)} rows, tables give {len(da)}")
for tr, r in zip(_trs, da):
    c = [re.sub(r"<[^>]+>", "", x).strip() for x in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)]
    c = c[-7:]
    ud = "—" if not r["n"] else (f"+{r['max_up']:.2f}" if r["max_up"] is not None else "—") + " / " + \
        (f"{MINUS}{r['max_down']:.2f}" if r["max_down"] is not None else "—")
    want = [r["label"], str(r["fitted"]), f"{r['n']} ({r['n_robust']})", ud,
            "—" if r["min_q"] is None else f"{r['min_q']:.1e}", num(r["median"]), r["chip"][0]]
    if c != want:
        misses.append(f"MaAsLin3 row {r['question'][:20]} / {r['label']}: report {c}, tables {want}")
    if r["n"] and not r["calls"].call_at_6.all():
        misses.append(f"{r['label']}: a call is not a call at >= 6 samples, but the heatmap caption says all are")
    for cl in r["calls"][r["calls"].robust].itertuples():
        expect(f"({num(cl.coef)}, q = {cl.q:.2g})", f"robust call {r['label']} {cl.taxon.split(';')[-1]}")

# ---- whole-design figures: taxon counts in the captions, images present ----
_taxa = sorted({t for r in da for t in r["calls"].taxon})
_gen = [t for t in _taxa if not species_named(t)]
expect(f"The {len(_taxa)} taxa with a call at q &lt; 0.05 in any row of the table, {len(_gen)} of them named only to genus or family", "heatmap taxa")
expect(f"limited to the {len(_gen)} taxa named no further than genus or family", "heatmap genus view")
expect("Every call is also a call at the stricter filter", "stricter filter claim")
for _f in ("factorial_heatmap", "factorial_forest", "gplus_factorial_heatmap", "gplus_factorial_forest"):
    if f'src="figures/MaAsLin3_factorial_Greengenes2_{_f}.png"' not in html:
        misses.append(f"report does not show {_f}")

# ---- interaction: joint q as the forest plot computes it, nominal count, scope ----
_ti = load_tables(TAB, "MaAsLin3_factorial_Greengenes2_int")
_qj = joint_q(_ti)
_rt = pd.read_csv(os.path.join(TAB, "maaslin3", "MaAsLin3_factorial_Greengenes2_min2_rows.tsv"), sep="\t")
_ri = _rt[(_rt.model == "abundance") & (_rt.family == "interaction") & _rt.p.notna()]
if abs(min(q for d in _qj.values() for q in d.values()) - _ri.q.min()) > 1e-12:
    misses.append("interaction: the forest plot's joint q and the rows table's q disagree")
if (_ri.q < 0.05).any():
    misses.append("interaction: a taxon has q < 0.05 but the text says none survives")
_nom, _tot = int((_ri.p < 0.05).sum()), len(_ri)
expect(f"there are {_nom} of them among {_tot} estimates, where chance alone would give about {round(_tot * 0.05)}", "nominal interaction count")
expect(f"{_nom} of {_tot} estimates nominal where about {round(_tot * 0.05)} are expected", "five-ways nominal count")
expect(f"has {_nom} nominal estimates where about {round(_tot * 0.05)} are", "finding nominal count")
expect(f"the smallest q is {_ri.q.min():.2g}", "smallest interaction q (caption)")
expect(f"(smallest q = {_ri.q.min():.2g})", "smallest interaction q (five ways, key finding)")
_per = _ri.groupby("row").size()
expect(f"{_per['Interaction: Inulin vs Cellulose']}, {_per['Interaction: RS vs Cellulose']} and {_per['Interaction: RS vs Inulin']} taxa", "forest scope")

# ---- presence model and quoted findings ----
_pv = _rt[(_rt.model == "prevalence") & _rt.p.notna()]
if (_pv.q < 0.05).any():
    misses.append("presence model: a call exists but the text says none")
expect(f"the smallest p is {_pv.p.min():.2g} and the smallest q {_pv.q.min():.2g}", "presence model floor")
R = {(r["question"], r["label"]): r for r in da}
d0 = lambda l: R[("Do the diets differ in unirradiated mice?", l)]
d1 = lambda l: R[("Do the diets differ in irradiated mice?", l)]
ir = lambda l: R[("Does irradiation change the microbiome within each diet?", l)]
_akk = d0("Inulin vs Cellulose")["calls"]; _akk = _akk[_akk.taxon.str.contains("Akkermansia")].iloc[0]
expect(f"up {_akk.coef:.2f} log2 (about a {round(2 ** _akk.coef, -1):.0f}-fold larger share)", "Akkermansia inulin")
expect(f"+{_akk.coef:.2f} for Inulin vs Cellulose means a taxon holding about a {round(2 ** _akk.coef, -1):.0f}-fold", "table note example")
expect(f"inulin changes the share of {d0('Inulin vs Cellulose')['n']} taxa against cellulose and RS changes {d0('RS vs Cellulose')['n']}", "finding diet counts")
expect(f"inulin and RS change the share of {d0('Inulin vs Cellulose')['n']} and {d0('RS vs Cellulose')['n']} taxa against cellulose, "
       f"{d0('Inulin vs Cellulose')['n_robust']} and {d0('RS vs Cellulose')['n_robust']} of them robustly", "synthesis diet counts")
expect(f"The typical taxon's share is {-d0('Inulin vs Cellulose')['median']:.2f} log2 lower on inulin", "inulin median")
expect(f"The two fibres differ from each other in {d0('RS vs Inulin')['n']} taxa, none robust", "fibre vs fibre")
_ar = d0("RS vs Inulin")["calls"]; _ar = _ar[_ar.taxon.str.contains("Akkermansia")].iloc[0]
expect(f"{-_ar.coef:.2f} log2 lower on RS than on inulin, which holds in {_ar.loo_holds_in} of the 19 refits", "Akkermansia RS vs inulin")
expect(f"Inulin+IR and RS+IR differ from Cellulose+IR in {d1('Inulin+IR vs Cellulose+IR')['n']} and {d1('RS+IR vs Cellulose+IR')['n']} taxa", "irradiated diet counts")
_i18 = int((d1("Inulin+IR vs Cellulose+IR")["calls"].loo_holds_in == 18).sum())
expect(f"five hold in 18 of the 19 refits" if _i18 == 5 else f"{_i18} hold in 18", "Inulin+IR near misses")
expect(f"({-d1('Inulin+IR vs Cellulose+IR')['median']:.2f} and {-d1('RS+IR vs Cellulose+IR')['median']:.2f} log2)", "irradiated medians")
expect(f"RS+IR against Inulin+IR gives {d1('RS+IR vs Inulin+IR')['n']} calls, none robust", "RS+IR vs Inulin+IR")
expect(f"({num(ir('Cellulose+IR vs Cellulose')['median'])} log2)", "cellulose irradiation median")
expect(f"the typical taxon's share falls {-ir('RS+IR vs RS')['median']:.2f} log2", "RS irradiation median")
expect(f"RS has the most calls ({ir('RS+IR vs RS')['n']})", "RS irradiation calls")
if any(ir(l)["n_robust"] for l in ("Cellulose+IR vs Cellulose", "Inulin+IR vs Inulin", "RS+IR vs RS")):
    misses.append("an irradiation row has a robust call, but the text says none")

# ---- structural zeros: shares of reads quoted in the text ----
import biom  # noqa: E402
from plot_untestable_heatmap import collapse  # noqa: E402

sh = pd.read_csv(os.path.join(TAB, f"{PRES}_untestable_shares.tsv"), sep="\t")
expect(f"{sh.untestable_pct.mean():.2f}% of an average sample's reads", "structural zeros mean share")
gm = sh.groupby("Group").untestable_pct.mean()
expect(f"{gm['RS+IR']:.1f}% in RS+IR samples, {gm['Cellulose']:.1f}% in Cellulose and {gm['Cellulose+IR']:.1f}% in Cellulose+IR", "structural zeros group shares")
rest = gm[["RS", "Inulin", "Inulin+IR"]]
expect(f"against {rest.min():.1f}–{rest.max():.1f}% in the RS and inulin groups", "structural zeros other-group shares")
if rest.max() >= gm[["RS+IR", "Cellulose", "Cellulose+IR"]].min():
    misses.append("structural zeros: the three named groups are not the three highest")
_z = untestable[[c for c in untestable.columns if " (n=" in c]].eq(0).sum(axis=1)
expect(f"{int((_z == 1).sum())} are missing from one group only, the other {int((_z > 1).sum())} from two or more", "structural zeros by number of groups")
fmeta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t", index_col=0, keep_default_na=False)
fc = biom.load_table(os.path.join(os.path.dirname(OUT), "tables", "study-seqs.biom")).to_dataframe(dense=True)[fmeta.index]
ftax = pd.read_csv(os.path.join(os.path.dirname(OUT), "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                   sep="\t", index_col=0)["Taxon"]
rel = collapse(fc, ftax) / fc.sum() * 100
row = lambda k: rel.loc[[t for t in rel.index if t.endswith(k)][0]]
am = row("Amulumruptor caecigallinarius"); g = fmeta.Group
if sh.loc[sh.top_taxon_pct.idxmax(), "top_taxon"] != am.name:
    misses.append("Amulumruptor is not the largest single share among the structural zeros")
ri = am[g == "RS+IR"]
expect(f"({ri.min():.0f}–{ri.max():.0f}%) of the reads in the three RS+IR samples, at most {am[g != 'RS+IR'].max():.2f}% anywhere else",
       "Amulumruptor shares")
if (ri == 0).any() or am[g.isin(["Inulin", "Inulin+IR"])].gt(0).any():
    misses.append("Amulumruptor presence does not match the text")
ms = row("Mucispirillum schaedleri")
if ms[g.isin(["Inulin", "Inulin+IR", "RS", "RS+IR"])].gt(0).any() or ms[g.isin(["Cellulose", "Cellulose+IR"])].eq(0).any():
    misses.append("Mucispirillum is not in every cellulose sample and no fibre sample")

# ---- structural zeros: presence table, carriers, complete separations, dot map ----
import check_structural_zeros  # noqa: E402

sz_miss, szdf = check_structural_zeros.check(html, text)
misses += sz_miss
misses += check_structural_zeros.check_dots(html, OUT)
for label, n_abs, n in (("RS+IR vs Cellulose+IR", 31, 36), ("RS+IR vs Inulin+IR", 22, 27), ("RS+IR vs RS", 19, 22)):
    d = szdf[szdf.contrast == label]
    if (int((d.absent_from == "RS+IR").sum()), len(d)) != (n_abs, n):
        misses.append(f"direction {label}: tables give {int((d.absent_from == 'RS+IR').sum())} of {len(d)}")
    expect(f"{n_abs} of the {n} against", f"direction {label}")
_cell = lambda k: [int(v) for v in untestable.set_index("taxon").loc[[t for t in untestable.taxon if t.endswith(k)][0]].drop("n_samples")]
for k, want in (("Mucispirillum schaedleri", [3, 0, 0, 4, 0, 0]), ("Bifidobacterium globosum", [0, 3, 2, 0, 4, 0]),
                ("g__Duncaniella;s__", [3, 3, 0, 4, 4, 0]), ("Lactobacillus johnsonii", [0, 3, 0, 4, 1, 0])):
    if _cell(k) != want:
        misses.append(f"{k}: presence per group {_cell(k)}, the text describes {want}")
# Structural zeros are no longer drawn per contrast; any cladogram left in the
# report must still match its sidecar.
fig_miss, _ = check_structural_zeros.check_figures(html, OUT)
misses += fig_miss
if 'untestable_heatmap.png"' in html:
    misses.append("the per-sample structural-zero heatmap is shown; the dot map replaces it")

# ---- taxa per sample by group (presence_interaction.R's by-cell table) ----
_bc = pd.read_csv(os.path.join(TAB, f"{PRES}_presence_by_cell.tsv"), sep="\t", keep_default_na=False)
_bc = {(r[0], r[1]): r.taxa_present for r in _bc.itertuples(index=False)}
for _d in ("Cellulose", "Inulin", "RS"):
    expect(f"{_bc[(_d, 'Control')]:.1f} to {_bc[(_d, 'IR')]:.1f}", f"taxa per sample {_d}")
expect(f"RS samples carry {_bc[('RS', 'Control')]:.1f} of the 84 taxa before irradiation and {_bc[('RS', 'IR')]:.1f}", "finding RS taxa per sample")

# ---- DA part order: opener, structural zeros, MaAsLin3, findings ----
_da = html.index("<h2>Differential abundance</h2>")
_order = [html.find(h, _da) for h in ("Structural zeros: taxa present", "Absent from first / second", "_presence_dots.png",
                                      "MaAsLin3: change in share of reads", "Typical taxon", "factorial_heatmap.png",
                                      "Interaction coefficients", "<h3>Findings</h3>", "<h2>Summary</h2>")]
if -1 in _order or _order != sorted(_order):
    misses.append(f"DA parts missing or out of order: {_order}")
for _gone in ("How taxa are tested", ">Both together<", "Reading the numbers"):
    if _gone in html:
        misses.append(f"'{_gone}' is still in the report")

# ---- community composition: group means of the species that differ most, shares only ----
_c = fc.loc[fc.index.intersection(ftax.index)]
_paths = ftax.loc[_c.index].map(lambda t: ";".join([p for p in (x.strip() for x in t.split(";")) if p][:7]))
_c = _c.groupby(_paths).sum()
_rel = _c / _c.sum() * 100
_G = ["Cellulose", "Inulin", "RS", "Cellulose+IR", "Inulin+IR", "RS+IR"]
_mean = pd.DataFrame({g: _rel.loc[:, fmeta.loc[_rel.columns, "Group"] == g].mean(axis=1) for g in _G})
_mean = _mean[~_mean.index.str.endswith("g__;s__")]
_top = (_mean.max(axis=1) - _mean.min(axis=1)).sort_values(ascending=False).index[:12]
_c0 = html.index("<h2>Community composition</h2>")
_cs = html[_c0:html.index("</section>", _c0)]
_rows = re.findall(r"<tr><td>.*?</tr>", _cs)
if len(_rows) != 12:
    misses.append(f"composition table has {len(_rows)} rows, expected 12")
_pct = lambda v: "&lt;0.1%" if v < 0.05 else f"{v:.1f}%"
for _r, _t in zip(_rows, _top):
    _got = re.findall(r'<td class="n">(.*?)</td>', _r)
    if _got != [_pct(v) for v in _mean.loc[_t, _G]]:
        misses.append(f"composition row {_t.split(';')[-1]}: report {_got}, tables {[_pct(v) for v in _mean.loc[_t, _G]]}")
_under = int((_mean.loc[_top, "Inulin"] < 0.1).sum())
expect(f"unirradiated inulin mice have almost none (under 0.1% of reads) of {_under} of the 12", "composition inulin absences")
expect(f"Unirradiated inulin mice have almost none of {_under} of the 12 species that differ most", "summary composition")
_sp = lambda k: _mean.loc[[t for t in _mean.index if t.endswith(k)][0]]
_bc_, _ak = _sp("caccae"), _sp("muciniphila_D_776786")
expect(f"({_bc_['Inulin']:.1f}% and {_bc_['Inulin+IR']:.1f}% of reads, unirradiated and irradiated)", "composition B. caccae")
expect(f"({_ak['Inulin']:.1f}% and {_ak['Inulin+IR']:.1f}%)", "composition Akkermansia")
_am = _sp("Amulumruptor caecigallinarius")
expect(f"{_am['RS+IR']:.1f}% of reads against at most {_am.drop('RS+IR').max():.2f}% in any other group", "composition Amulumruptor")
_ub = _sp("UBA7173 sp002491305")
expect(f"sp002491305 ({_ub['RS+IR']:.1f}%) with unirradiated RS ({_ub['RS']:.1f}%)", "composition UBA7173")
_es = _sp("g__Escherichia")
expect(f"({_es['Cellulose']:.1f}% and {_es['Cellulose+IR']:.1f}%)", "composition Escherichia")
if _es.idxmax() not in ("Cellulose",) or _es["Cellulose+IR"] < _es[["Inulin", "Inulin+IR", "RS+IR"]].max():
    misses.append("composition: cellulose groups do not carry the most Escherichia")
_mu = _sp("Mucispirillum schaedleri")
expect(f"({_mu['Cellulose']:.1f}% unirradiated)", "composition Mucispirillum")
if (_mu.drop(["Cellulose", "Cellulose+IR"]) > 0).any():
    misses.append("composition: Mucispirillum is not cellulose-only")
_ph = fc.loc[fc.index.intersection(ftax.index)]
_ph = _ph.groupby(ftax.loc[_ph.index].map(lambda t: t.split(";")[1].strip() if ";" in t else "unassigned")).sum()
_ph = _ph / _ph.sum() * 100
_bd = pd.Series({g: _ph.loc["p__Bacteroidota", fmeta.index[fmeta.Group == g]].mean() for g in _G})
expect(f"rises from {_bd['Cellulose']:.1f}% of reads on cellulose to {_bd['RS+IR']:.1f}% in irradiated RS mice", "composition Bacteroidota")
if _bd.idxmin() != "Cellulose" or _bd.idxmax() != "RS+IR":
    misses.append("composition: Bacteroidota is not lowest on cellulose and highest in RS+IR")
if re.search(r"\bq\b|\bq =|raw p", re.sub(r"<[^>]+>", " ", _cs)):
    misses.append("Community composition quotes a test result; it precedes the DA methods")
for _rk in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"):
    _v = f"rk-{_rk[0].lower()}"
    if (f'<div class="view {_v}">' not in _cs or f"taxa_factorial_Greengenes2_{_rk}_groups.png" not in _cs
            or f'value="{_v}"' not in _cs or f".view.{_v}{{display:block}}" not in html):
        misses.append(f"composition rank panel: {_rk} view missing or not wired")
if 'value="rk-g" checked' not in _cs:
    misses.append("composition rank panel: Genus is not the default view")
_order2 = [html.find(h) for h in ("<h2>Beta diversity</h2>", "<h2>Community composition</h2>", "<th>Species · Greengenes2</th>",
                                  "The twelve species whose", '<figure class="panel ranks">', "per-sample bars", "<h2>Differential abundance</h2>")]
if -1 in _order2 or _order2 != sorted(_order2):
    misses.append(f"composition parts missing or out of order: {_order2}")

# ---- Predicted function (PICRUSt2; run_picrust2.sh, picrust2_permanova.R, run_picrust2_maaslin3.R) ----
import html as _html  # noqa: E402
_P0 = html.index("<h2>Predicted function (PICRUSt2)</h2>")
_ps = html[_P0:html.index("</section>", _P0)]
_pt = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", _ps))
_PIPE = os.path.join(OUT, "picrust2", "pipeline")
_FILES = {"MetaCyc": "pathways_out/path_abun_unstrat.tsv.gz", "EC": "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
          "KO": "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"}
_ab = {lv: pd.read_csv(os.path.join(_PIPE, f), sep="\t", index_col=0)[fmeta.index] for lv, f in _FILES.items()}
_npred = {lv: int((a.sum(axis=1) > 0).sum()) for lv, a in _ab.items()}
expect(f"{_npred['MetaCyc']:,} MetaCyc pathways, {_npred['EC']:,} enzymes (EC numbers) and {_npred['KO']:,} KEGG orthologs (KOs) are predicted", "PICRUSt2 feature counts")
_nsti = pd.read_csv(os.path.join(_PIPE, "combined_marker_predicted_and_nsti.tsv.gz"), sep="\t", index_col=0).metadata_NSTI
_wn = pd.read_csv(os.path.join(_PIPE, "KO_metagenome_out", "weighted_nsti.tsv.gz"), sep="\t", index_col=0).weighted_NSTI[fmeta.index]
if (_nsti > 2).any():
    misses.append("PICRUSt2: an ASV exceeds NSTI 2, but the text says none is dropped")
_gm = _wn.groupby(fmeta.Group).median()
expect(f"ran on the study's {len(_nsti):,} ASVs", "PICRUSt2 ASV count")
expect(f"(median {_nsti.median():.3f}, largest {_nsti.max():.2f})", "NSTI median and max")
expect(f"runs from {_wn.min():.3f} to {_wn.max():.3f} per sample", "weighted NSTI range")
expect(f"medians run from {_gm.min():.3f} ({_gm.idxmin()}) to {_gm.max():.3f} ({_gm.idxmax()})", "weighted NSTI group medians")
_pm = pd.read_csv(os.path.join(TAB, "PICRUSt2_permanova.tsv"), sep="\t").set_index(["level", "term"])
_pp = lambda p: "0.0001" if p <= 1e-4 else (f"{p:.4f}" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}"))
for _lv, _name in (("MetaCyc", "MetaCyc pathways"), ("EC", "Enzymes (EC)"), ("KO", "KEGG orthologs (KO)")):
    _r = lambda t: _pm.loc[(_lv, t)]
    expect(f"{_name} " + " ".join(f"{_r(t).effect * 100:.1f}% ({_pp(_r(t).p)})" for t in ("Diet", "Condition", "Diet:Condition"))
           + f" {_pp(_r('group').p)}", f"PICRUSt2 PERMANOVA row {_lv}")
_rng = lambda t, f: (min(f(_pm.loc[(lv, t)]) for lv in _FILES), max(f(_pm.loc[(lv, t)]) for lv in _FILES))
for _t, _lab in (("Diet", "Diet explains"), ("Condition", "Irradiation explains")):
    lo, hi = _rng(_t, lambda r: r.effect * 100)
    expect(f"{_lab} {lo:.0f}–{hi:.0f}%", f"PICRUSt2 {_t} R2 range")
_plo, _phi = _rng("Condition", lambda r: r.p)
expect(f"(p = {_pp(_plo)}–{_pp(_phi)})", "PICRUSt2 irradiation p range")
_ilo, _ihi = _rng("Diet:Condition", lambda r: r.effect * 100)
_iplo, _iphi = _rng("Diet:Condition", lambda r: r.p)
expect(f"diet-by-irradiation term {_ilo:.0f}–{_ihi:.0f}% (p = {_pp(_iplo)}–{_pp(_iphi)})", "PICRUSt2 interaction")
if _iplo < 0.05:
    misses.append("PICRUSt2 interaction PERMANOVA now passes, but the text says it stays short")
expect(f"PERMDISP p ≥ {_pp(min(_pm.loc[(lv, 'group')].p for lv in _FILES))}", "PICRUSt2 PERMDISP")
# The PCA outlier named in the captions: the sample farthest out on PC1 of the scaled PCA, every database.
for _lv, _a in _ab.items():
    _a = _a[_a.std(axis=1) > 0]
    _z = ((_a.T - _a.T.mean()) / _a.T.std()).values
    _u, _sv, _ = np.linalg.svd(_z - _z.mean(0), full_matrices=False)
    _pc = pd.Series(np.abs(_u[:, 0] * _sv[0]), index=fmeta.index)
    if _pc.idxmax() != "d26T1.10" or fmeta.Group[_pc.idxmax()] != "Cellulose":
        misses.append(f"PICRUSt2 PCA {_lv}: the far-out sample is {_pc.idxmax()}, the caption says d26T1.10")
if _pt.count("One unirradiated cellulose sample, d26T1.10, sits far out on the first axis in all three databases") != 3:
    misses.append("PICRUSt2 PCA captions do not each name d26T1.10 as the far-out sample")
# MaAsLin3 table, every cell, from the rows tables.
_RD = [("Inulin vs Cellulose", "Control"), ("RS vs Cellulose", "Control"), ("RS vs Inulin", "Control"),
       ("Inulin vs Cellulose", "IR"), ("RS vs Cellulose", "IR"), ("RS vs Inulin", "IR"),
       ("IR within Cellulose", ""), ("IR within Inulin", ""), ("IR within RS", ""),
       ("Interaction: Inulin vs Cellulose", ""), ("Interaction: RS vs Cellulose", ""), ("Interaction: RS vs Inulin", "")]
_MS = {}
for _k, _lv in (("mc", "MetaCyc"), ("ec", "EC"), ("ko", "KO")):
    _d = pd.read_csv(os.path.join(TAB, f"PICRUSt2_{_lv}_maaslin3_rows.tsv"), sep="\t", keep_default_na=False, na_values=["NA"])
    _a = _d[(_d.model == "abundance") & _d.p.notna()]
    _v0 = _ps.index(f'<div class="view {_k}">', _ps.index("MaAsLin3 calls per row"))
    _vt = _ps[_v0:_ps.index("</table>", _v0)]
    _trs = re.findall(r"<tr[^>]*>(.*?)</tr>", _vt[_vt.index("<tbody>"):], re.S)
    if len(_trs) != 12:
        misses.append(f"PICRUSt2 {_lv} table: {len(_trs)} rows, expected 12")
    for _tr, (_row, _st) in zip(_trs, _RD):
        _x = _a[(_a.row == _row) & (_a.stratum == _st)]
        _rb = _x[_x.robust]
        _MS[(_lv, _row, _st)] = dict(robust=len(_rb), median=_x.coef.median(), down=int((_rb.coef < 0).sum()))
        _want = [f"{len(_x):,}", f"{int(_x.call.sum()):,}", f"{len(_rb):,}",
                 "—" if not len(_rb) else f"{int((_rb.coef > 0).sum())} / {int((_rb.coef < 0).sum())}", num(_x.coef.median())]
        _got = [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", _tr, re.S)][-5:]
        if _got != _want:
            misses.append(f"PICRUSt2 {_lv} {_row} {_st}: report {_got}, tables {_want}")
    _pv = _d[(_d.model == "prevalence") & _d.q.notna()]
    if _pv.call.any():
        misses.append(f"PICRUSt2 {_lv}: the presence model has a call, but the text says none")
_M = lambda lv, r, st="": _MS[(lv, r, st)]
expect(f"RS+IR differs from Cellulose+IR in {_M('MetaCyc', 'RS vs Cellulose', 'IR')['robust']} pathways, {_M('EC', 'RS vs Cellulose', 'IR')['robust']:,} enzymes and {_M('KO', 'RS vs Cellulose', 'IR')['robust']:,} KOs, where the typical pathway's share is {-_M('MetaCyc', 'RS vs Cellulose', 'IR')['median']:.2f} log2 lower", "PICRUSt2 RS+IR vs Cellulose+IR")
expect(f"inulin differs from cellulose in {_M('KO', 'Inulin vs Cellulose', 'Control')['robust']} KOs and RS in {_M('KO', 'RS vs Cellulose', 'Control')['robust']}", "PICRUSt2 unirradiated diets")
expect(f"Irradiation changes {_M('MetaCyc', 'IR within RS')['robust']} pathways, {_M('EC', 'IR within RS')['robust']} enzymes and {_M('KO', 'IR within RS')['robust']} KOs on RS", "PICRUSt2 irradiation on RS")
expect(f"({_M('KO', 'IR within RS')['down']} of the {_M('KO', 'IR within RS')['robust']} KOs)", "PICRUSt2 RS direction")
expect(f"but only {_M('MetaCyc', 'IR within Cellulose')['robust']} pathways on cellulose and none on inulin", "PICRUSt2 cellulose and inulin")
if any(_M(lv, "IR within Inulin")["robust"] for lv in _FILES):
    misses.append("PICRUSt2: irradiation within inulin has a robust call, but the text says none")
_irs, _iri = _M("MetaCyc", "Interaction: RS vs Cellulose"), _M("MetaCyc", "Interaction: RS vs Inulin")
expect(f"The interaction rows hold {_irs['robust']} and {_iri['robust']} robust pathways for RS against cellulose and inulin, all lower on RS", "PICRUSt2 interaction counts")
if _irs["down"] != _irs["robust"] or _iri["down"] != _iri["robust"]:
    misses.append("PICRUSt2: an RS interaction call is higher on RS, but the text says all are lower")
_mc = pd.read_csv(os.path.join(TAB, "PICRUSt2_MetaCyc_maaslin3_rows.tsv"), sep="\t", keep_default_na=False, na_values=["NA"])
_mc = _mc[(_mc.model == "abundance") & _mc.robust]
_mc["description"] = _mc.description.map(_html.unescape)
_ferm = ["mixed acid fermentation", "superpathway of anaerobic sucrose degradation", "glycolysis III (from glucose)"]
_fr = [r for r in ("IR within RS", "Interaction: RS vs Cellulose", "Interaction: RS vs Inulin")]
_fv = _mc[_mc.row.isin(_fr) & _mc.description.isin(_ferm)]
if len(_fv) != 9 or (_fv.coef >= 0).any():
    misses.append(f"PICRUSt2: the three fermentation pathways are not robust and lower in all three RS rows ({len(_fv)} of 9)")
_fmed = [_M("MetaCyc", r)["median"] for r in _fr]
expect(f"{-_fv.coef.max():.2f} to {-_fv.coef.min():.2f} log2 lower, against a typical shift of {num(max(_fmed))} to {num(min(_fmed))} log2", "PICRUSt2 fermentation sizes")
_ser = _mc[(_mc.row == "IR within Cellulose")]
if sorted(_ser.description) != ["L-serine biosynthesis II", "formaldehyde assimilation I (serine pathway)"] or (_ser.coef <= 0).any():
    misses.append(f"PICRUSt2: the cellulose irradiation calls are {sorted(_ser.description)}")
for _r in _ser.itertuples():
    expect(f"{_r.description} ({num(_r.coef)}", f"PICRUSt2 serine {_r.description}")
_rel = _ab["MetaCyc"] / _ab["MetaCyc"].sum() * 100
for _id in _ser.id:
    if _rel.loc[_id, fmeta.Group == "Cellulose+IR"].min() <= _rel.loc[_id, fmeta.Group == "Cellulose"].max():
        misses.append(f"PICRUSt2: {_id} is not higher in every irradiated cellulose sample, but the text says so")
for _lv, _k in (("MetaCyc", "mc"), ("EC", "ec"), ("KO", "ko")):
    for _f in ("pca", "heatmap", "volcano"):
        if f'src="figures/picrust2_{_lv}_{_f}.png"' not in _ps:
            misses.append(f"PICRUSt2: {_lv} {_f} figure not shown")
_pord = [_ps.find(x) for x in ("PICRUSt2 predicts which genes", "Prediction coverage</h3>", "picrust2_nsti.png", "Whole-profile tests</h3>",
                               "Diet: R² (p)", "Diet explains", "picrust2_MetaCyc_pca.png", "Per-feature tests</h3>", "<strong>MaAsLin3.</strong>",
                               "MaAsLin3 calls per row", "Robust calls gather", "picrust2_MetaCyc_heatmap.png", "picrust2_MetaCyc_volcano.png",
                               "Two limits apply", "Short-chain fatty acid and formate metabolism</h3>", "picrust2_acetate.png",
                               "picrust2_butyrate.png", "Butyrate synthesis moves most", "<h3>Findings</h3>")]
if -1 in _pord or _pord != sorted(_pord):
    misses.append(f"PICRUSt2 section out of order or missing a part: {_pord}")
if not (html.index("<h2>Differential abundance</h2>") < _P0 < html.index("<h2>Summary</h2>")):
    misses.append("PICRUSt2 section is not between Differential abundance and Summary")
if _ps.count('<div class="finding') != 4:
    misses.append("PICRUSt2 section should close on four findings")
if re.search(r"LinDA|Fisher|Firth", _pt):
    misses.append("PICRUSt2 section mentions a method it no longer uses")

# ---- PICRUSt2 PCA ellipses: plot_ggpicrust2.R draws them for every group with >= 3 samples ----
_small = sorted(meta.Group.value_counts().loc[lambda n: n < 3].index)
if _small != ["RS"]:
    misses.append(f"PCA captions say every group but unirradiated RS has an ellipse; groups under 3 samples are {_small}")
if html.count("Dashed ellipses span one standard deviation of each group with at least three samples (every group but unirradiated RS)") != 3:
    misses.append("PCA captions: the ellipse sentence is not in all three")

# ---- PICRUSt2 SCFA and formate panels: every quoted range from plot_picrust2_targets.py's tables ----
_SC = {_m: pd.read_csv(os.path.join(TAB, f"PICRUSt2_{_m}_results.tsv"), sep="\t", dtype={"feature": str})
       for _m in ("acetate", "formate", "propionate", "butyrate")}
for _m in _SC:
    if f'src="figures/picrust2_{_m}.png"' not in html:
        misses.append(f"PICRUSt2: {_m} panel not shown")


def _sg(v, nd):
    return (MINUS if v < 0 else "+") + f"{abs(v):.{nd}f}"


def _rg(vals, nd=1):
    vals = sorted(vals, key=abs)
    a, b = _sg(vals[0], nd), _sg(vals[-1], nd)
    return a if a == b else f"{a} to {b}"


def _robust(m, rows, ids, n, sign=-1, why=""):
    d = _SC[m][_SC[m].row.isin(rows) & _SC[m].feature.isin(ids)]
    if len(d) != n or not d.robust.all() or not (d.coef * sign > 0).all():
        misses.append(f"PICRUSt2 {m}: {why} are not {n} robust calls in the stated direction")
    return d.coef


def _typ(m, rows, dbs):
    return _SC[m][_SC[m].row.isin(rows) & _SC[m].database.isin(dbs)].typical_shift.unique()


_INU, _RSC, _RSI = ["Inulin vs Cellulose", "Inulin+IR vs Cellulose+IR"], ["RS+IR vs Cellulose+IR"], ["RS+IR vs Inulin+IR"]
_inu = _rg(_robust("butyrate", _INU, ["2.3.1.9", "4.2.1.17", "ko:K00626", "ko:K01715"], 8, why="inulin thiolase/crotonase"))
_inu_ty = _rg(_typ("butyrate", _INU, ["EC", "KO"]))
_four = _rg(_robust("butyrate", _RSC, ["2.3.1.9", "1.1.1.157", "4.2.1.17", "1.3.8.1", "ko:K00626", "ko:K00074",
                                       "ko:K01715", "ko:K00248"], 8, why="RS+IR butyrate enzymes"))
_rsc_ty = _rg(_typ("butyrate", _RSC, ["EC", "KO"]))
expect(f"thiolase and crotonase are robust calls at {_inu} log2, where the typical feature shifts {_inu_ty}", "SCFA inulin butyrate")
expect(f"against Cellulose+IR all four are robust calls at {_four} log2, against a typical shift of {_rsc_ty}", "SCFA RS+IR butyrate")
_glut = _rg(_robust("butyrate", _RSC, ["2.8.3.12", "ko:K01039", "ko:K01040", "ko:K01615"], 4, why="glutarate route"))
_kal = _sg(_robust("butyrate", _RSC, ["ko:K18014"], 1, why="kal").iloc[0], 1)
expect(f"the glutarate route ({_glut}) and the lysine route's last enzyme ({_kal})", "SCFA glutarate and lysine")
if _SC["butyrate"][_SC["butyrate"].feature.isin(["ko:K00634", "ko:K00929", "2.3.1.19", "2.7.2.7"])].robust.any():
    misses.append("PICRUSt2 butyrate: ptb/buk has a robust call, but the text says none")
expect("butyrate kinase with phosphate butyryltransferase, has no robust call in any row", "SCFA butyrate final step")
_pfl_i = _rg(_robust("formate", _RSI, ["2.3.1.54", "ko:K00656", "ko:K04069"], 3, why="pfl vs Inulin+IR"), 2)
_mix_i = _sg(_robust("formate", _RSI, ["FERMENTATION-PWY"], 1, why="mixed acid fermentation").iloc[0], 2)
expect(f"activating enzyme pflA are robust calls at {_pfl_i} log2", "SCFA pfl vs Inulin+IR")
expect(f"at {_mix_i}, where the typical feature shifts {_rg(_typ('formate', _RSI, ['MetaCyc', 'EC', 'KO']))}", "SCFA mixed acid vs Inulin+IR")
_pfl_c = _rg(_robust("formate", _RSC, ["2.3.1.54", "1.97.1.4", "ko:K00656", "ko:K04069", "ko:K04070"], 5, why="pfl vs Cellulose+IR"))
expect(f"features are robust calls at {_pfl_c}, but the typical feature there shifts {_rg(_typ('formate', _RSC, ['MetaCyc', 'EC', 'KO']))}",
       "SCFA pfl vs Cellulose+IR")
expect(f"is a robust call at {_rg(_robust('formate', _RSC + _RSI, ['6.3.4.3', 'ko:K01938'], 4, why='fhs'))} against both other irradiated diets",
       "SCFA fhs")
expect(f"(interaction {_rg(_robust('formate', ['Interaction: RS vs Inulin'], ['6.3.4.3', 'ko:K01938'], 2, why='fhs interaction'))})",
       "SCFA fhs interaction")
_mut = _rg(_robust("propionate", _RSC, ["5.4.99.2", "5.1.99.1", "ko:K01847", "ko:K05606"], 4, sign=1, why="mutase/epimerase"), 2)
expect(f"methylmalonyl-CoA mutase and epimerase are robust calls at {_mut} log2 against a typical shift of {_rsc_ty}", "SCFA propionate mutase")
expect(f"methylmalonyl-CoA carboxyltransferase are at {_rg(_robust('propionate', _RSC, ['P108-PWY', '2.1.3.1', 'ko:K03416'], 3, why='succinate pathway'))}",
       "SCFA propionate pathway")
_pa = _SC["acetate"][(_SC["acetate"].role == "Acetate from acetyl-CoA (Pta–AckA)") & _SC["acetate"].robust]
if len(_pa) != 2 or _pa.coef.abs().max() >= 0.4:
    misses.append(f"PICRUSt2 acetate: Pta–AckA has {len(_pa)} robust calls, largest {_pa.coef.abs().max():.2f}; the text says 2, under 0.4")
expect(f"has {len(_pa)} robust calls, both under 0.4 log2", "SCFA Pta-AckA")
_ato = _rg(_robust("acetate", ["Inulin+IR vs Cellulose+IR"] + _RSC, ["ko:K01034", "ko:K01035"], 4, why="AtoDA"))
expect(f"is {_ato} in Inulin+IR and RS+IR against Cellulose+IR", "SCFA AtoDA")
expect(f"thiolase and crotonase are robust calls at {_inu} log2, unirradiated and irradiated", "SCFA finding inulin")
expect(f"all four enzymes are robust calls at {_four}, against {_rsc_ty}", "SCFA finding RS+IR")

# ---- every section: at most four findings, and they close it ----
# The Summary alone may add up to two red findings ("finding no": easily misread results) after its
# green or orange ones, with only a heading between the two groups (user, 2026-10-01).
for _m in re.finditer(r"<section[^>]*>(.*?)</section>", html, re.S):
    _sec = _m.group(1)
    _h = re.search(r"<h2>(.*?)</h2>", _sec)
    _name = _h.group(1) if _h else "(section without a heading)"
    _boxes = list(re.finditer(r'<div class="finding([^"]*)">.*?</div>', _sec, re.S))
    _red = [b for b in _boxes if " no" in b.group(1)] if _name == "Summary" else []
    _other = [b for b in _boxes if b not in _red]
    if len(_other) > 4:
        misses.append(f"{_name}: {len(_other)} findings, at most 4")
    if len(_red) > 2:
        misses.append(f"{_name}: {len(_red)} red findings, at most 2")
    if _red and _other and _red[0].start() < _other[-1].start():
        misses.append(f"{_name}: a red finding comes before a green/orange one")
    if _boxes:
        _after = re.sub(r"\s+", " ", re.sub(r"<script.*?</script>|<[^>]+>", " ", _sec[_boxes[-1].end():], flags=re.S)).strip()
        if _after:
            misses.append(f"{_name}: content after its findings: {_after[:60]!r}")
        _rest = re.sub(r'<script.*?</script>|<div class="finding[^"]*">.*?</div>', "", _sec[_boxes[0].start():], flags=re.S)
        if _red and _other:   # a heading may introduce the red group
            _rest = re.sub(r"<h3[^>]*>.*?</h3>", "", _rest, flags=re.S)
        _between = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", _rest)).strip()
        if _between:
            misses.append(f"{_name}: findings are not together at the end ({_between[:60]!r} between them)")
    if _name == "Differential abundance" and len(_boxes) < 3:
        misses.append(f"Differential abundance: {len(_boxes)} findings, 3 to 4 close the section")

# LinDA is not used in this report and is not mentioned (user, 2026-10-04).
if "LinDA" in html:
    misses.append("the report mentions LinDA, which it does not use")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {n_fig} figures, all quoted numbers found")
