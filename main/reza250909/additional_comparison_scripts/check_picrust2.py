"""Report check for reza250909's Predicted function (PICRUSt2) section, MaAsLin3 version (2026-10-04).

check(html) recomputes, from picrust2/ and the tables, the feature counts, NSTI, every PERMANOVA
cell, every MaAsLin3 table cell, the counts and directions the prose quotes, the curated SCFA
numbers, the Summary row and the section's order, independently of the script that wrote them.
Returns (misses, ko_nominal, ko_chance) for the key-findings check.
"""
import html as _html
import os
import re

import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
TAB = os.path.join(OUT, "tables")
PIPE = os.path.join(OUT, "picrust2", "pipeline")
FILES = {"MetaCyc": "pathways_out/path_abun_unstrat.tsv.gz", "KO": "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
         "EC": "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"}
ROWS = ["SCFAs vs Control", "SCFAs+IR vs IR", "IR vs Control", "SCFAs+IR vs SCFAs", "Interaction: SCFAs vs Control"]
MINUS = "−"
num = lambda v: f"+{v:.2f}" if v >= 0 else f"{MINUS}{-v:.2f}"
pp = lambda p: "0.0001" if p <= 1e-4 else (f"{p:.4f}" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}"))
B = "SCFAs+IR vs IR"


def check(html):
    misses = []
    p0 = html.index("<h2>Predicted function (PICRUSt2)</h2>")
    sec = html[p0:html.index("</section>", p0)]
    text = _html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", sec)))

    def expect(v, why, where=None):
        if v not in (where or text):
            misses.append(f"PICRUSt2: {why}: '{v}' not in report")

    meta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t", index_col=0)
    ab = {lv: pd.read_csv(os.path.join(PIPE, f), sep="\t", index_col=0)[meta.index] for lv, f in FILES.items()}
    npred = {lv: int((a.sum(axis=1) > 0).sum()) for lv, a in ab.items()}
    expect(f"{npred['KO']:,} KEGG orthologs (KO), {npred['EC']:,} enzymes (EC) and {npred['MetaCyc']:,} MetaCyc pathways", "feature counts")
    nsti = pd.read_csv(os.path.join(PIPE, "combined_marker_predicted_and_nsti.tsv.gz"), sep="\t", index_col=0).metadata_NSTI
    wn = pd.read_csv(os.path.join(PIPE, "KO_metagenome_out", "weighted_nsti.tsv.gz"), sep="\t", index_col=0).weighted_NSTI[meta.index]
    if (nsti > 2).any():
        misses.append("PICRUSt2: an ASV exceeds NSTI 2, but the text says none is dropped")
    expect(f"{len(nsti):,} ASVs in 32 mice", "ASV count")
    expect(f"(median {nsti.median():.3f}, largest {nsti.max():.2f})", "NSTI")
    expect(f"weighted NSTI runs from {wn.min():.2f} to {wn.max():.2f}", "weighted NSTI")

    pm = pd.read_csv(os.path.join(TAB, "PICRUSt2_permanova.tsv"), sep="\t").set_index(["level", "term"])
    for lv, name in (("MetaCyc", "MetaCyc pathways"), ("KO", "KEGG orthologs (KO)"), ("EC", "Enzymes (EC)")):
        expect(f"{name} " + " ".join(f"{pm.loc[(lv, t)].effect * 100:.1f}% ({pp(pm.loc[(lv, t)].p)})"
                                     for t in ("Diet", "Radiation", "Diet:Radiation")) + f" {pp(pm.loc[(lv, 'group')].p)}",
               f"PERMANOVA row {lv}")
    rng = lambda t: f"{min(pm.loc[(lv, t)].effect for lv in FILES) * 100:.1f}–{max(pm.loc[(lv, t)].effect for lv in FILES) * 100:.1f}%"
    prng = lambda t: f"p = {pp(min(pm.loc[(lv, t)].p for lv in FILES))}–{pp(max(pm.loc[(lv, t)].p for lv in FILES))}"
    expect(f"diet explains {rng('Diet')} of the variation ({prng('Diet')})", "diet range")
    expect(f"the interaction {rng('Diet:Radiation')} ({prng('Diet:Radiation')})", "interaction range")
    if min(pm.loc[(lv, t)].p for lv in FILES for t in ("Diet", "Radiation", "Diet:Radiation")) < 0.05:
        misses.append("PICRUSt2: a whole-profile term passes, but the text says none does")

    S = {}
    calls = []
    for k, lv in (("mc", "MetaCyc"), ("ko", "KO"), ("ec", "EC")):
        d = pd.read_csv(os.path.join(TAB, f"PICRUSt2_{lv}_maaslin3_rows.tsv"), sep="\t", keep_default_na=False, na_values=["NA"])
        a = d[(d.model == "abundance") & d.p.notna()]
        v0 = sec.index(f'<div class="view {k}">', sec.index("MaAsLin3 per row"))
        vt = sec[v0:sec.index("</table>", v0)]
        trs = re.findall(r"<tr[^>]*>(.*?)</tr>", vt[vt.index("<tbody>"):], re.S)
        if len(trs) != len(ROWS):
            misses.append(f"PICRUSt2 {lv} table: {len(trs)} rows, expected {len(ROWS)}")
        for tr, row in zip(trs, ROWS):
            x = a[a.row == row]
            nom = x[x.p < 0.05]
            S[(lv, row)] = dict(n=len(x), nom=len(nom), ch=round(0.05 * len(x)), lo=int((nom.coef < 0).sum()),
                                hi=int((nom.coef > 0).sum()), med=x.coef.median(), minq=x.q.min())
            calls += [(lv, row, r) for r in x[x.call].itertuples()]
            want = [f"{len(x):,}", f"{len(nom):,} ({round(0.05 * len(x)):,})", f"{int((nom.coef < 0).sum())} / {int((nom.coef > 0).sum())}",
                    f"{int(x.call.sum())} ({int(x.robust.sum())})", num(x.coef.median())]
            got = [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)][-5:]
            if got != want:
                misses.append(f"PICRUSt2 {lv} {row}: report {got}, tables {want}")
        if d[(d.model == "prevalence") & d.q.notna()].call.any():
            misses.append(f"PICRUSt2 {lv}: the presence model has a call, but the text says none")
    s = lambda lv, row: S[(lv, row)]
    expect(f"in {s('MetaCyc', B)['nom']} pathways, {s('KO', B)['nom']:,} KOs and {s('EC', B)['nom']} enzymes at raw p < 0.05, "
           f"where chance gives about {s('MetaCyc', B)['ch']}, {s('KO', B)['ch']} and {s('EC', B)['ch']}, in both directions "
           f"({s('KO', B)['lo']} KOs lower, {s('KO', B)['hi']} higher)", "SCFAs+IR vs IR excess")
    if not all(s(lv, B)["nom"] > 2.5 * s(lv, B)["ch"] for lv in FILES) or any(
            s(lv, r)["nom"] >= 1.5 * s(lv, r)["ch"] for lv in FILES for r in ROWS if r != B):
        misses.append("PICRUSt2: the excess is no longer confined to SCFAs+IR vs IR")
    expect(f"(smallest q {s('KO', B)['minq']:.2f} among KOs)", "SCFAs+IR vs IR smallest q")
    if len(calls) != 1 or calls[0][2].robust:
        misses.append(f"PICRUSt2: {len(calls)} calls at q < 0.05 (robust: {[c[2].robust for c in calls]}); the text says one, not robust")
    else:
        lv, row, c = calls[0]
        expect(f"{_html.unescape(c.description)} in {row} ({num(c.coef)} log2, q = {c.q:.2f})", "the one call")
        expect(f"({int(c.loo_holds_in)} of 32)", "the one call's refits")
    expect(f"within {max(abs(v['med']) for v in S.values()):.2f} log2 of zero", "typical shift bound")

    # Curated SCFA panels, recomputed from the MaAsLin3 rows (not the panels' own tables).
    ko = pd.read_csv(os.path.join(TAB, "PICRUSt2_KO_maaslin3_rows.tsv"), sep="\t", keep_default_na=False, na_values=["NA"])
    ko = ko[ko.model == "abundance"].set_index(["id", "row"])
    but = [ko.loc[(f, B)] for f in ("ko:K00626", "ko:K00074", "ko:K01715", "ko:K00248", "ko:K00634", "ko:K00929")]
    if not all(b.p < 0.05 and b.coef < 0 for b in but):
        misses.append("PICRUSt2: a butyrate-synthesis KO is not lower at raw p < 0.05 in SCFAs+IR vs IR")
    expect(f"({num(max(b.coef for b in but))} to {num(min(b.coef for b in but))} log2, raw p < 0.05), where the typical KO shifts {num(s('KO', B)['med'])}",
           "butyrate synthesis")
    fct, oxc = ko.loc[("ko:K07749", B)], ko.loc[("ko:K01577", B)]
    if not (fct.p < 0.05 and fct.coef < 0 and oxc.p < 0.05 and oxc.coef > 0):
        misses.append("PICRUSt2: formyl-CoA transferase / oxalyl-CoA decarboxylase no longer move apart")
    expect(f"formyl-CoA transferase is lower ({num(fct.coef)}) and oxalyl-CoA decarboxylase higher ({num(oxc.coef)})", "formate")
    pr = pd.read_csv(os.path.join(TAB, "PICRUSt2_propionate_results.tsv"), sep="\t")
    if (pr.p < 0.05).any():
        misses.append("PICRUSt2: a propionate feature reaches raw p < 0.05, but the text says none")
    for m in ("acetate", "formate", "propionate", "butyrate"):
        if f'src="figures/picrust2_{m}.png"' not in sec:
            misses.append(f"PICRUSt2: {m} panel not shown")
    for lv in FILES:
        for w in ("pca", "heatmap", "volcano"):
            if f'src="figures/picrust2_{lv}_{w}.png"' not in sec:
                misses.append(f"PICRUSt2: {lv} {w} not shown")
    order = [sec.find(x) for x in ("PICRUSt2 predicts which genes", "Prediction coverage</h3>", "picrust2_nsti.png",
                                   "Whole-profile tests</h3>", "Diet: R² (p)", "In pathways, orthologs and enzymes alike",
                                   "picrust2_MetaCyc_pca.png", "Per-feature tests</h3>", "<strong>MaAsLin3.</strong>",
                                   "MaAsLin3 per row", "One row stands out", "picrust2_MetaCyc_heatmap.png",
                                   "picrust2_MetaCyc_volcano.png", "Short-chain fatty acid and formate metabolism</h3>",
                                   "picrust2_acetate.png", "Butyrate synthesis is where", "No curated feature passes correction",
                                   "What the three databases show together</h3>", 'class="finding')]
    if -1 in order or order != sorted(order):
        misses.append(f"PICRUSt2 section out of order or missing a part: {order}")
    if sec.count('<div class="finding') != 4:
        misses.append("PICRUSt2 section should close on four findings")
    summ = _html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html[html.index("<h2>Summary</h2>"):])))
    expect(f"differ from irradiated control-diet mice in {s('KO', B)['nom']:,} KOs at raw p < 0.05 against about {s('KO', B)['ch']} by chance",
           "Summary row", summ)
    i = s("KO", "Interaction: SCFAs vs Control")
    expect(f"{s('MetaCyc', 'Interaction: SCFAs vs Control')['nom']} pathways, {i['nom']} KOs and {s('EC', 'Interaction: SCFAs vs Control')['nom']} enzymes at raw p < 0.05, "
           f"below the {s('MetaCyc', 'Interaction: SCFAs vs Control')['ch']}, {i['ch']} and {s('EC', 'Interaction: SCFAs vs Control')['ch']} chance gives", "interaction row", summ)
    if not all(s(lv, "Interaction: SCFAs vs Control")["nom"] < s(lv, "Interaction: SCFAs vs Control")["ch"] for lv in FILES):
        misses.append("PICRUSt2: the per-feature interaction is no longer below chance in every database")
    return misses, s("KO", B)["nom"], s("KO", B)["ch"]
