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


def load(name):
    return pd.read_csv(os.path.join(TAB, name), sep="\t").set_index("taxon")


def hits(df, coef="TimepointDay25"):
    ok = (df[f"q_{coef}"] < 0.05) & df[f"passed_ss_{coef}"].astype(str).str.upper().eq("TRUE")
    return df[ok]


per_arm = {a: load(f"ANCOMBC2_day17day25_{a}_Greengenes2.tsv") for a in ("IFA", "IFA_IR")}
unpaired = {a: load(f"ANCOMBC2_day17day25_{a}_Greengenes2_unpaired.tsv") for a in ("IFA", "IFA_IR")}
joint = {a: load(f"ANCOMBC2_day17day25_joint_Greengenes2_armref{a}.tsv") for a in ("IFA", "IFA_IR")}

expect(f"{per_arm['IFA'].q_TimepointDay25.min():.3f} / {joint['IFA'].q_TimepointDay25.min():.4f}", "IFA min q pair")
expect(f"{per_arm['IFA_IR'].q_TimepointDay25.min():.3f} / {joint['IFA_IR'].q_TimepointDay25.min():.4f}", "IFA+IR min q pair")
est = {a: d.se_TimepointDay25.notna().sum() for a, d in per_arm.items()}
expect(f"{len(per_arm['IFA'])} taxa in IFA and {len(per_arm['IFA_IR'])} in IFA+IR", "per-arm tested")
expect(f"Of those, {est['IFA']} and {est['IFA_IR']} have an estimable effect", "per-arm estimable")
expect(f"return {len(hits(unpaired['IFA']))} taxa in IFA and {len(hits(unpaired['IFA_IR']))} in IFA+IR", "unpaired counts")
both = hits(per_arm["IFA_IR"]).index.intersection(hits(joint["IFA_IR"]).index)
# DA table counts: per-arm, joint, both.
for arm, label in (("IFA", "IFA"), ("IFA_IR", "IFA+IR")):
    cells = (len(hits(per_arm[arm])), len(hits(joint[arm])),
             len(hits(per_arm[arm]).index.intersection(hits(joint[arm]).index)))
    want = "".join(f'<td class="n">{c}</td>' for c in cells)
    if f"<td>{label}</td>{want}" not in html:
        misses.append(f"DA table row {label}: tables give per-arm/joint/both {cells}")
fail = lambda d: int(((d.q_TimepointDay25 < 0.05) & ~d.index.isin(hits(d).index)).sum())
expect(f"The per-arm fit finds no taxon: {fail(per_arm['IFA_IR'])} reach q &lt; 0.05", "IFA+IR per-arm failing count")
expect(f"{fail(joint['IFA_IR'])} more reach q &lt; 0.05 there", "IFA+IR joint failing count")
expect(f"has {fail(joint['IFA'])} taxa at q &lt; 0.05 with IFA as the reference arm, and all {fail(joint['IFA'])} fail",
       "IFA joint failing count")
if len(hits(per_arm["IFA_IR"])) or len(hits(per_arm["IFA"])) or len(hits(joint["IFA"])) or len(both):
    misses.append("DA prose says no per-arm hit, no joint IFA hit and no taxon in both models")
hj = hits(joint["IFA_IR"]).lfc_TimepointDay25 / LN2
names = {t.split(";")[-2 if t.rstrip().endswith("__") else -1].strip(): v for t, v in hj.items()}
if sorted(names) != ["g__Borkfalkia", "s__Lawsonibacter sp000492175"]:
    misses.append(f"joint IFA+IR hits are {sorted(names)}, prose names Borkfalkia and Lawsonibacter")
else:
    expect(f"+{names['g__Borkfalkia']:.2f} log2 (about {2 ** names['g__Borkfalkia']:.0f}-fold", "Borkfalkia estimate")
    expect(f"sp000492175, +{names['s__Lawsonibacter sp000492175']:.2f}", "Lawsonibacter estimate")
expect(f"tests {len(joint['IFA'])} taxa", "joint tested")
ic = "TimepointDay25:ArmIFA_IR"
expect(f"smallest q is {joint['IFA'][f'q_{ic}'].min():.3f}", "interaction min q")
sing = pd.read_csv(os.path.join(TAB, "ANCOMBC2_day17day25_joint_Greengenes2_armrefIFA_singular.tsv"), sep="\t")
expect(f"{int(sing.singular.sum())} of the joint model's {len(sing)}", "joint singular count")

# Composition table: mean relative abundance per arm and timepoint, summed to
# ANCOM-BC2's species-level paths, next to each arm's per-arm q.
import biom  # noqa: E402

STUDY = os.path.dirname(THIS_DIR)
counts = biom.load_table(os.path.join(STUDY, "tables", "study-seqs.biom")).to_dataframe(dense=True)
tax = pd.read_csv(os.path.join(STUDY, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                  sep="\t", index_col=0)["Taxon"]
counts = counts.loc[counts.index.intersection(tax.index)]
paths = tax.loc[counts.index].map(lambda t: ";".join([p for p in (x.strip() for x in t.split(";")) if p][:7]))
counts = counts.groupby(paths).sum()
ra = {}
for a in ("IFA", "IFA_IR"):
    meta = pd.read_csv(os.path.join(TAB, f"metadata_subset_{a}.tsv"), sep="\t", index_col=0)
    rel = counts[meta.index] / counts[meta.index].sum()
    for d in ("Day17", "Day25"):
        ra[a, d] = rel[meta.index[meta.Timepoint == d]].mean(axis=1) * 100
delta = sum(ra[a, "Day25"] - ra[a, "Day17"] for a in ("IFA", "IFA_IR")) / 2
delta = delta[~delta.index.str.endswith("g__;s__")]
top = delta.abs().sort_values(ascending=False).index[:11]


def pct(v):
    return "0%" if v == 0 else ("&lt;0.1%" if v < 0.05 else f"{v:.1f}%")


def qcell(a, taxon):
    if taxon not in per_arm[a].index or pd.isna(per_arm[a].q_TimepointDay25.get(taxon)):
        return "—"
    q = per_arm[a].q_TimepointDay25[taxon]
    return f"{q:.3f}" + ("*" if q < 0.05 and taxon not in hits(per_arm[a]).index else "")


for taxon in top:
    row = " ".join(c for a in ("IFA", "IFA_IR")
                   for c in (pct(ra[a, "Day17"][taxon]), pct(ra[a, "Day25"][taxon]), qcell(a, taxon)))
    expect(row, f"composition row {taxon.split(';')[-1]}")
expect(f"In IFA+IR, {len(hits(per_arm['IFA_IR']).index.intersection(top))} of the 11 are per-arm hits: "
       f"{sum(per_arm['IFA_IR'].q_TimepointDay25.get(t, 1) < 0.05 for t in top)} reach q", "composition IFA+IR hit count")
if len(hits(per_arm["IFA"]).index.intersection(top)):
    misses.append("composition: prose says no IFA per-arm hits among the 11")
ms = per_arm["IFA"].index[per_arm["IFA"].index.str.contains("schaedleri")][0]
expect(f"joint q {joint['IFA'].q_TimepointDay25[ms]:.3f} in IFA, {joint['IFA_IR'].q_TimepointDay25[ms]:.3f}",
       "M. schaedleri joint q")
for phylum in ("p__Campylobacterota_A", "p__Deferribacterota"):
    present = counts[(counts.index.str.split(";").str[1] == phylum)].sum(axis=1).gt(0).sum()
    if present != 1:
        misses.append(f"composition: {phylum} has {present} species-level taxa, prose says one")

# ---- untestable-taxa heatmap caption ----
import biom  # noqa: E402
from plot_untestable_heatmap import collapse  # noqa: E402

stem = "ANCOMBC2_day17day25_joint_Greengenes2"
un = pd.read_csv(os.path.join(TAB, f"{stem}_untestable_absent_from_a_cell.tsv"), sep="\t")
sh = pd.read_csv(os.path.join(TAB, f"{stem}_untestable_shares.tsv"), sep="\t")
cc = {c.split(" (n=")[0]: c for c in un.columns if "(n=" in c}
absent = {k: un[c].eq(0) for k, c in cc.items()}
only = lambda *ks: sum(absent[k] for k in cc).eq(len(ks)) & np.logical_and.reduce([absent[k] for k in ks])
n_ifa, n_ir, n_both = (int(only(*k).sum()) for k in (("IFA Day17",), ("IFA+IR Day17",), ("IFA Day17", "IFA+IR Day17")))
d17_only = n_ifa + n_ir + n_both
expect(f"The {len(un)} taxa that pass the presence rule", "heatmap taxon count")
expect(f"{d17_only} of the {len(un)} are missing only at Day 17, {n_ifa} of them in the IFA arm alone, "
       f"{n_ir} in IFA+IR alone and {n_both} in both arms", "Day-17-only counts")
gm = sh.groupby("Group").untestable_pct.mean()
others = gm.drop("IFA Day17")
expect(f"{gm['IFA Day17']:.2f}% of reads on average, against {others.min():.1f}–{others.max():.1f}%", "heatmap cell shares")
fmeta = pd.read_csv(os.path.join(OUT, "metadata_joint_cells.tsv"), sep="\t", index_col=0, keep_default_na=False)
fc = biom.load_table(os.path.join(os.path.dirname(OUT), "tables", "study-seqs.biom")).to_dataframe(dense=True)[fmeta.index]
ftax = pd.read_csv(os.path.join(os.path.dirname(OUT), "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                   sep="\t", index_col=0)["Taxon"]
rel = collapse(fc, ftax) / fc.sum() * 100
row = lambda k: rel.loc[[t for t in rel.index if t.endswith(k)][0]]
g = fmeta.Group
al = row("Alloprevotella sp002933955")
if al[g == "IFA Day17"].gt(0).any() or al[g != "IFA Day17"].eq(0).any():
    misses.append("heatmap: Alloprevotella is not present in exactly the 12 non-IFA-Day17 samples")
expect(f"in all 12 other samples ({al[g != 'IFA Day17'].min():.1f}–{al.max():.1f}%)", "Alloprevotella range")
ac = row("Bacteroides_H_857956 acidifaciens")
if ac[g.str.startswith("IFA Day")].gt(0).any():
    misses.append("heatmap: B. acidifaciens is seen in the IFA arm")
d17, d25 = ac[g == "IFA+IR Day17"], ac[g == "IFA+IR Day25"]
expect(f"rises from {d17.min():.0f}–{d17.max():.1f}% of reads at Day 17 to {d25.min():.1f}–{d25.max():.1f}% at Day 25",
       "B. acidifaciens ranges")
ut = set(un.taxon)
for arm, f, word in (("IFA", "ANCOMBC2_day17day25_IFA_Greengenes2.tsv", "IFA"),
                     ("IFA_IR", "ANCOMBC2_day17day25_IFA_IR_Greengenes2.tsv", "IFA+IR")):
    pa = pd.read_csv(os.path.join(TAB, f), sep="\t")
    pa = pa[pa.taxon.isin(ut)]
    se = [c for c in pa if c.startswith("se_") and "Intercept" not in c][0]
    q = "q_" + se[3:]
    est = pa[pa[se].notna()]
    hit = est[(est[q] < 0.05) & est["passed_ss_" + se[3:]].astype(str).str.upper().eq("TRUE")]
    expect(f"{len(est)} of these taxa in IFA" if word == "IFA" else f"and {len(est)} in IFA+IR", f"per-arm SE count {word}")
    if len(hit):
        misses.append(f"heatmap: {len(hit)} per-arm {word} hits among untestable taxa, caption says none")

# ---- structural zeros: presence table and finding ----
import check_structural_zeros  # noqa: E402

check_structural_zeros.sz.configure("day17-day25")
sz_miss, szdf = check_structural_zeros.check(html, text)
misses += sz_miss
floor = check_structural_zeros.sz.mcnemar_exact(4, 0)
expect(f"cannot go below p = {floor:.3f}", "paired floor")
for label, grp, n_abs, n in (("IFA", "IFA Day17", 31, 35), ("IFA+IR", "IFA+IR Day17", 16, 22)):
    d = szdf[szdf.contrast == label]
    if (int((d.absent_from == grp).sum()), len(d)) != (n_abs, n):
        misses.append(f"direction {label}: tables give {int((d.absent_from == grp).sum())} of {len(d)}")
    expect(f"{n_abs} of the {n}", f"direction {label}")

fig_miss, n_views = check_structural_zeros.check_figures(html, OUT)
misses += fig_miss
if n_views == 0:
    misses.append("no view draws structural zeros")

# ---- presence turnover (presence_turnover_day17_day25.R) ----
_s = pd.read_csv(os.path.join(TAB, "ANCOMBC2_day17day25_joint_Greengenes2_presence_turnover.tsv"), sep="\t").iloc[0]
_per = pd.read_csv(os.path.join(TAB, "ANCOMBC2_day17day25_joint_Greengenes2_presence_turnover_by_mouse.tsv"), sep="\t")
expect(f"Mean turnover is {_s.turnover_a:.2f} in IFA and {_s.turnover_b:.2f} in IFA+IR (exact p = {_s.p_exact:.2f})", "turnover and p")
expect(f"Turnover {_s.turnover_a:.2f} vs {_s.turnover_b:.2f}, exact p = {_s.p_exact:.2f}", "turnover table row")
expect(f"IFA mice gain {_s.gained_a:.1f} of these taxa on average and lose {_s.lost_a:.1f}; IFA+IR mice gain {_s.gained_b:.1f} and lose {_s.lost_b:.1f}", "gained and lost")
_r = lambda a: f"{_per.day17[_per.arm == a].min()}–{_per.day17[_per.arm == a].max()}"
expect(f"with {_r('IFA')} of the {_s.taxa} taxa at Day 17 against {_r('IFA+IR')} in IFA+IR", "Day 17 ranges")
expect(f"is asked of all {_s.taxa} taxa at once: the {_s.tested} the joint model tests and the {_s.untestable} it leaves out", "turnover taxon set")
expect(f"{_s.taxa} taxa pass the rule in at least one arm: {_s.tested} are seen in each of the four", "joint taxon split")
expect(f"The plot covers the {_s.tested} taxa seen in all four arm × timepoint cells. The {_s.untestable} absent", "forest scope")
# The text says these are the factorial report's unirradiated IFA mice at Day 17.
_fm = pd.read_csv(os.path.join(os.path.dirname(OUT), "metadata.tsv"), sep="\t", index_col=0, keep_default_na=False)
_fifa = set(_fm.index[(_fm.Timepoint == "Day17") & (_fm.Treatment == "IFA") & (_fm.Radiation == "None")])
_jm = pd.read_csv(os.path.join(OUT, "metadata_joint_cells.tsv"), sep="\t", index_col=0)
if _fifa != set(_jm.index[_jm.Group == "IFA Day17"]):
    misses.append("the IFA Day 17 samples are not the factorial report's unirradiated IFA group")

# ---- summary section ----
_fl = beta.permanova_q.max()
if not (beta.permanova_q == beta.permanova_q.min()).all():
    misses.append("summary says every PERMANOVA row sits at the permutation floor; they differ")
expect(f"(R² {beta.R2.min() * 100:.0f}–{beta.R2.max() * 100:.0f}%)", "summary R2 range")
_spread = beta[beta.permdisp_q < 0.05].groupby("arm").metric.apply(set)
if not all(v == {"jaccard", "braycurtis", "weighted_unifrac"} for v in _spread):
    misses.append(f"summary names Jaccard, Bray-Curtis and weighted UniFrac as the spread metrics; tables give {dict(_spread)}")
expect(f"in IFA+IR {fail(per_arm['IFA_IR'])} reach q &lt; 0.05 but fail the sensitivity check", "summary per-arm failing count")
expect(f"The joint model finds {len(hits(joint['IFA_IR']))} in IFA+IR", "summary joint count")
expect(f"(0 of {len(joint['IFA'])} taxa) or in presence (turnover p = {_s.p_exact:.2f})", "summary arm comparison")

# ---- structural-zero evidence (structural_zero_evidence.py) ----
_S = os.path.join(TAB, "ANCOMBC2_day17day25_joint_Greengenes2")
_cal = pd.read_csv(_S + "_structural_zeros_calibration.tsv", sep="\t").set_index("contrast")
_dep = pd.read_csv(_S + "_structural_zeros_depth.tsv", sep="\t")
_t = _cal.loc["all contrasts"]
expect(f"{_t.observed} calls in {_dep.taxon.nunique()} taxa show complete separation between the timepoints, where "
       f"swapping timepoints within mice gives {_t.expected:.2f} on average (exact p = {_t.p_perm:.3f})", "SZ calibration")
expect(f"cannot go below p = {_cal.loc['IFA'].p_perm:.3f}", "SZ per-arm floor")
if (_dep.q_depth >= 0.05).any():
    misses.append("text says no complete-separation absence is explained by depth")
_app = _dep.present_in.str.endswith("Day25").sum()
expect(f"{['No', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight'][_app]} of the {['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten'][len(_dep)]} calls are taxa that appear", "SZ appear count")
expect(f"{_t.observed} calls against {_t.expected:.2f} expected by chance", "SZ summary numbers")
_rows = re.findall(r"<tr><td>(.*?)</td><td>(IFA|IFA\+IR)</td><td>(appears by Day 25|gone by Day 25)</td><td class=\"n\">([0-9.]+)%</td></tr>", html)
_want = sorted((r.contrast, "appears by Day 25" if r.present_in.endswith("Day25") else "gone by Day 25", f"{r.mean_share_pct:.2f}")
               for r in _dep.itertuples())
if sorted((r[1], r[2], r[3]) for r in _rows) != _want:
    misses.append("SZ table rows do not match the tables")
_al = _dep[_dep.taxon.str.endswith("Alloprevotella sp002933955")]
if list(zip(_al.contrast, _al.present_in)) != [("IFA", "IFA Day25")]:
    misses.append("Alloprevotella is not the IFA arm's appear-by-Day-25 call")
_de = _dep[_dep.taxon.str.endswith("Desulfovibrio sp011039135")]
if list(zip(_de.contrast, _de.present_in)) != [("IFA+IR", "IFA+IR Day17")]:
    misses.append("Desulfovibrio is not the IFA+IR arm's lost-by-Day-25 call")

# ---- restructured DA section: synthesis claims ----
expect(f"mice gain {_s.gained_a:.1f} (IFA) and {_s.gained_b:.1f} (IFA+IR) of the {_s.taxa} taxa on average", "synthesis gains")
expect(f"({_r('IFA')} at Day 17 against {_r('IFA+IR')} in IFA+IR)", "synthesis Day 17 ranges")
expect(f"{_t.observed} calls against {_t.expected:.2f} expected by chance), and mice gain", "synthesis calibration")
_un = pd.read_csv(os.path.join(TAB, "ANCOMBC2_day17day25_joint_Greengenes2_untestable_absent_from_a_cell.tsv"), sep="\t")
_ab = _un[[c for c in _un.columns if "(n=" in c]].eq(0).rename(columns=lambda c: c.split(" (n=")[0])
_d17 = _ab["IFA Day17"] & ~_ab["IFA Day25"] & ~_ab["IFA+IR Day25"] & ~_ab["IFA+IR Day17"]
expect(f"{int(_d17.sum())} of the 33 taxa missing only at Day 17 are missing in the IFA arm alone", "synthesis IFA-only count")
expect(f"0 of {len(joint['IFA'])} taxa, smallest q {joint['IFA'][f'q_{ic}'].min():.3f}", "synthesis interaction")
for _h in ("How taxa are tested", "Structural zeros: taxa that appear or vanish", "ANCOM-BC2: abundance in the taxa seen", "Both together"):
    if _h not in html:
        misses.append(f"DA part heading missing: {_h}")
_order = [html.find(_h) for _h in ("How taxa are tested", "Structural zeros: taxa that appear", "ANCOM-BC2: abundance in the taxa", "Both together", "<h2>Summary</h2>")]
if _order != sorted(_order):
    misses.append("DA parts are out of order")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {len(re.findall(r'src=', html))} figures, all quoted numbers found")
