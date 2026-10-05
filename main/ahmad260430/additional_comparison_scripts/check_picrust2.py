"""Report check for the Predicted function (PICRUSt2) section of both ahmad260430 reports.

check(html, report) recomputes, from picrust2/ and the tables, the feature counts, NSTI, every
PERMANOVA cell, every MaAsLin3 table cell, the targeted-panel numbers the prose quotes and the
section's order, independently of the script that wrote the section. Returns a list of misses.
Called by check_report_factorial.py and check_report_day17_day25.py.
"""
import html as _html
import os
import re

import numpy as np
import pandas as pd

from picrust2_targets import FEATURES

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(THIS_DIR)
PIPE = os.path.join(STUDY, "picrust2", "pipeline")
FILES = {"MetaCyc": "pathways_out/path_abun_unstrat.tsv.gz", "EC": "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
         "KO": "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"}
MINUS = "−"
num = lambda v: f"+{v:.2f}" if v >= 0 else f"{MINUS}{-v:.2f}"
pp = lambda p: "0.0001" if p <= 1e-4 else (f"{p:.4f}" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}"))
REPORTS = {
    "factorial": dict(dir="comparison-Treatment-x-Radiation-factorial", rows=[
        ("IFA vs Control", "None", "IFA vs Control"), ("Inosine vs Control", "None", "Inosine vs Control"),
        ("IFA vs Inosine", "None", "IFA vs Inosine"), ("IFA vs Control", "IR", "IFA+IR vs IR"),
        ("Inosine vs Control", "IR", "Inosine+IR vs IR"), ("IFA vs Inosine", "IR", "IFA+IR vs Inosine+IR"),
        ("IR within Control", "", "IR vs Control"), ("IR within IFA", "", "IFA+IR vs IFA"),
        ("IR within Inosine", "", "Inosine+IR vs Inosine"), ("Interaction: IFA vs Control", "", "Interaction: IFA vs Control"),
        ("Interaction: Inosine vs Control", "", "Interaction: Inosine vs Control"),
        ("Interaction: IFA vs Inosine", "", "Interaction: IFA vs Inosine")]),
    "day17-day25": dict(dir="comparison-day17-day25", rows=[
        ("Day 17 to 25 in IFA", "", "Day 17 to 25, IFA"), ("Day 17 to 25 in IFA+IR", "", "Day 17 to 25, IFA+IR"),
        ("Average Day 17 to 25 shift", "", "Average Day 17 to 25 shift"),
        ("IFA+IR vs IFA at Day 17", "", "IFA+IR vs IFA, Day 17"), ("IFA+IR vs IFA at Day 25", "", "IFA+IR vs IFA, Day 25"),
        ("Average IFA+IR vs IFA", "", "Average IFA+IR vs IFA"),
        ("Interaction (IFA+IR shift minus IFA shift)", "", "Interaction (shift difference)")]),
}
AEROBIC = [(db, fid) for db, fid, role, _l in FEATURES["isoferulic"]
           if role.startswith(("CoA activation", "Ring")) or fid in ("ko:K21802", "1.2.1.67", "ko:K03862", "ko:K03863", "1.14.13.82")]


def check(html, report):
    misses = []
    cfg = REPORTS[report]
    out = os.path.join(STUDY, cfg["dir"])
    tab = os.path.join(out, "tables")
    p0 = html.index("<h2>Predicted function (PICRUSt2)</h2>")
    sec = html[p0:html.index("</section>", p0)]
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", sec))

    def expect(v, why):
        if v not in text:
            misses.append(f"PICRUSt2 {report}: {why}: '{v}' not in report")

    if report == "factorial":
        meta = pd.read_csv(os.path.join(tab, "metadata_day17.tsv"), sep="\t", index_col=0)
        samples = list(meta.index)
        group = meta.Group
    else:
        meta = pd.read_csv(os.path.join(out, "metadata_joint_cells.tsv"), sep="\t", index_col=0)
        samples = list(meta.index)
        group = meta.Group
    ab = {lv: pd.read_csv(os.path.join(PIPE, f), sep="\t", index_col=0)[samples] for lv, f in FILES.items()}
    npred = {lv: int((a.sum(axis=1) > 0).sum()) for lv, a in ab.items()}
    expect(f"{npred['MetaCyc']:,} MetaCyc pathways, {npred['EC']:,} enzymes (EC numbers) and {npred['KO']:,} KEGG orthologs (KOs) are predicted",
           "feature counts")
    nsti = pd.read_csv(os.path.join(PIPE, "combined_marker_predicted_and_nsti.tsv.gz"), sep="\t", index_col=0).metadata_NSTI
    wn = pd.read_csv(os.path.join(PIPE, "KO_metagenome_out", "weighted_nsti.tsv.gz"), sep="\t", index_col=0).weighted_NSTI[samples]
    if (nsti > 2).any():
        misses.append(f"PICRUSt2 {report}: an ASV exceeds NSTI 2, but the text says none is dropped")
    expect(f"(median {nsti.median():.3f}, largest {nsti.max():.2f})", "NSTI")
    expect(f"runs from {wn.min():.3f} to {wn.max():.3f} per sample", "weighted NSTI")

    pm = pd.read_csv(os.path.join(tab, "PICRUSt2_permanova.tsv"), sep="\t")
    if report == "factorial":
        pmi = pm.set_index(["level", "term"])
        for lv, name in (("MetaCyc", "MetaCyc pathways"), ("EC", "Enzymes (EC)"), ("KO", "KEGG orthologs (KO)")):
            expect(f"{name} " + " ".join(f"{pmi.loc[(lv, t)].effect * 100:.1f}% ({pp(pmi.loc[(lv, t)].p)})"
                                         for t in ("Treatment", "Radiation", "Treatment:Radiation"))
                   + f" {pp(pmi.loc[(lv, 'group')].p)}", f"PERMANOVA row {lv}")
        for t, lab in (("Treatment", "Treatment explains"), ("Radiation", "radiation")):
            e = [pmi.loc[(lv, t)].effect * 100 for lv in FILES]
            p = [pmi.loc[(lv, t)].p for lv in FILES]
            expect(f"{lab} {min(e):.0f}–{max(e):.0f}%" + (" of the variation in predicted function" if t == "Treatment" else "")
                   + f" (p = {pp(min(p))}–{pp(max(p))})", f"{t} range")
        if min(pmi.loc[(lv, "Radiation")].p for lv in FILES) < 0.05:
            misses.append("PICRUSt2 factorial: radiation PERMANOVA passes, but a finding says below threshold")
    else:
        pmi = pm.set_index(["level", "arm", "analysis"])
        for lv, name in (("MetaCyc", "MetaCyc pathways"), ("EC", "Enzymes (EC)"), ("KO", "KEGG orthologs (KO)")):
            for i, a in enumerate(("IFA", "IFA+IR")):
                r, d = pmi.loc[(lv, a, "permanova_bray_curtis")], pmi.loc[(lv, a, "permdisp_bray_curtis")]
                expect(("" if i else f"{name} ") + f"{a} {r.effect * 100:.1f}% ({pp(r.p)}) {pp(d.p)}", f"PERMANOVA {lv} {a}")
        e = pm[pm.analysis == "permanova_bray_curtis"].effect * 100
        expect(f"the day explains {e.min():.0f}–{e.max():.0f}% of the variation", "day R2 range")
        if (pm[(pm.analysis == "permdisp_bray_curtis") & (pm.arm == "IFA")].p > 1e-4).any():
            misses.append("PICRUSt2 day17-day25: IFA PERMDISP is not 0.0001 in every database, as the text says")

    # MaAsLin3 table cells, every database.
    res = {}
    for k, lv in (("mc", "MetaCyc"), ("ec", "EC"), ("ko", "KO")):
        d = pd.read_csv(os.path.join(tab, f"PICRUSt2_{lv}_maaslin3_rows.tsv"), sep="\t", keep_default_na=False, na_values=["NA"])
        a = d[(d.model == "abundance") & d.p.notna()]
        v0 = sec.index(f'<div class="view {k}">', sec.index("MaAsLin3 calls per row"))
        vt = sec[v0:sec.index("</table>", v0)]
        trs = re.findall(r"<tr[^>]*>(.*?)</tr>", vt[vt.index("<tbody>"):], re.S)
        got = {re.sub(r"<[^>]+>", "", c[-6]).strip(): c[-5:] for c in
               ([x for x in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)] for tr in trs)}
        got = {lab: [re.sub(r"<[^>]+>", "", x).strip() for x in v] for lab, v in got.items()}
        if len(trs) != len(cfg["rows"]):
            misses.append(f"PICRUSt2 {report} {lv} table: {len(trs)} rows, expected {len(cfg['rows'])}")
        for row, st, lab in cfg["rows"]:
            x = a[(a.row == row) & (a.stratum == st)]
            rb = x[x.robust]
            res[(lv, lab)] = dict(robust=len(rb), calls=int(x.call.sum()), median=x.coef.median())
            want = [f"{len(x):,}", f"{int(x.call.sum()):,}", f"{len(rb):,}",
                    "—" if not len(rb) else f"{int((rb.coef > 0).sum())} / {int((rb.coef < 0).sum())}", num(x.coef.median())]
            if got.get(lab) != want:
                misses.append(f"PICRUSt2 {report} {lv} {lab}: report {got.get(lab)}, tables {want}")
        pv = d[(d.model == "prevalence") & d.q.notna()]
        if pv.call.any():
            misses.append(f"PICRUSt2 {report} {lv}: the presence model has a call, but the text says none")
    R = lambda lv, lab: res[(lv, lab)]

    # Targeted panels, recomputed from the MaAsLin3 rows (not the panel's own results table).
    ko = pd.read_csv(os.path.join(tab, "PICRUSt2_KO_maaslin3_rows.tsv"), sep="\t", keep_default_na=False, na_values=["NA"])
    ec = pd.read_csv(os.path.join(tab, "PICRUSt2_EC_maaslin3_rows.tsv"), sep="\t", keep_default_na=False, na_values=["NA"])
    def feat(tbl, fid, row, st=""):
        x = tbl[(tbl.model == "abundance") & (tbl.id == fid) & (tbl.row == row) & (tbl.stratum == st)]
        return x.iloc[0]
    carr = lambda lv, fid, ss: int((ab[lv].loc[fid, ss] > 0).sum()) if fid in ab[lv].index else 0
    if report == "factorial":
        expect(f"IFA changes {R('MetaCyc', 'IFA vs Control')['robust']} pathways, {R('EC', 'IFA vs Control')['robust']:,} enzymes and "
               f"{R('KO', 'IFA vs Control')['robust']:,} KOs, where the typical pathway's share is {-R('MetaCyc', 'IFA vs Control')['median']:.2f} log2 lower",
               "IFA per-feature counts")
        brk = [feat(ko, f, "IFA vs Control", "None") for f in ("ko:K03784", "ko:K01239", "ko:K01488", "ko:K00087")]
        if not all(b.robust and b.coef < R("KO", "IFA vs Control")["median"] - 1.5 for b in brk):
            misses.append("PICRUSt2 factorial: an inosine-breakdown gene is not robust and well below the typical KO under IFA")
        c = [b.coef for b in brk]
        expect(f"robustly {-max(c):.2f} to {-min(c):.2f} log2 lower, where the typical KO is {-R('KO', 'IFA vs Control')['median']:.2f} log2 lower",
               "inosine breakdown under IFA")
        nuc = [feat(ko, f, "IFA vs Control", "None") for f in ("ko:K01081", "ko:K03787", "ko:K11751")]
        if not all(n.robust and n.coef > 0 for n in nuc):
            misses.append("PICRUSt2 factorial: a 5'-nucleotidase is not robustly higher under IFA")
        expect(f"{num(min(n.coef for n in nuc))} to {num(max(n.coef for n in nuc))} log2", "5'-nucleotidases")
        ino = ko[(ko.model == "abundance") & (ko.row == "Inosine vs Control") & (ko.stratum == "None")]
        curated = {f for db, f, _r, _l in FEATURES["inosine"] if db == "KO"}
        if ino[ino.id.isin(curated)].call.any():
            misses.append("PICRUSt2 factorial: Inosine changes a curated inosine KO, but the text says none")
        ir = [feat(ko, f, "Inosine vs Control", "IR") for f in ("ko:K03784", "ko:K01239", "ko:K01488")]
        if not all(x.robust and x.coef < 0 for x in ir):
            misses.append("PICRUSt2 factorial: Inosine+IR vs IR does not lower the three breakdown genes robustly")
        expect(f"({num(max(x.coef for x in ir))} to {num(min(x.coef for x in ir))} log2, typical KO {num(R('KO', 'Inosine+IR vs IR')['median'])})",
               "Inosine+IR breakdown genes")
        aer = max(carr(db, f, samples) for db, f in AEROBIC)
        expect(f"ring-cleavage enzymes are predicted in {aer} of the 24 mice, and phenolic acid decarboxylase in {carr('KO', 'ko:K13727', samples)}",
               "aerobic carriers")
        ccr = feat(ec, "1.3.1.108", "IFA vs Control", "None")
        if not ccr.robust or carr("EC", "1.3.1.108", samples) != len(samples):
            misses.append("PICRUSt2 factorial: caffeoyl-CoA reductase is not robust, or not predicted in every mouse")
        expect(f"is {-ccr.coef:.2f} log2 lower under IFA than in control mice and {-feat(ec, '1.3.1.108', 'IFA vs Control', 'IR').coef:.2f} lower with irradiation",
               "caffeoyl-CoA reductase")
    else:
        avg = "Average Day 17 to 25 shift"
        expect(f"Averaged over both arms, {R('MetaCyc', avg)['robust']} pathways, {R('EC', avg)['robust']} enzymes and {R('KO', avg)['robust']:,} KOs change",
               "averaged shift counts")
        for lab in ("IFA+IR vs IFA, Day 17", "IFA+IR vs IFA, Day 25", "Average IFA+IR vs IFA", "Interaction (shift difference)"):
            if any(R(lv, lab)["robust"] for lv in FILES):
                misses.append(f"PICRUSt2 day17-day25: {lab} has a robust call, but the text says none")
        up = [feat(ec, "3.2.2.1", avg)] + [feat(ko, f, avg) for f in ("ko:K03784", "ko:K09913", "ko:K03317")]
        if not all(x.robust and x.coef > R("KO", avg)["median"] for x in up):
            misses.append("PICRUSt2 day17-day25: an inosine-breakdown gene is not robustly above the typical shift")
        expect(f"robustly {num(min(x.coef for x in up))} to {num(max(x.coef for x in up))} log2 higher, where the typical KO rises {num(R('KO', avg)['median'])} log2",
               "inosine genes by Day 25")
        d17 = [s for s in samples if group[s].endswith("Day17")]
        d25 = [s for s in samples if group[s].endswith("Day25")]
        key = [("KO", "ko:K12508"), ("KO", "ko:K03862"), ("KO", "ko:K00448")]
        on25 = {s for db, f in key for s in d25 if carr(db, f, [s])}
        if any(carr(db, f, d17) for db, f in key):
            misses.append("PICRUSt2 day17-day25: an aerobic route is predicted at Day 17, but the text says none")
        n_ifa = sum(group[s].startswith("IFA Day") for s in on25)
        expect(f"predicted in {len(on25)} of the eight Day 25 samples ({n_ifa} IFA, {len(on25) - n_ifa} IFA+IR)", "Day 25 aerobic carriers")
        iso = pd.read_csv(os.path.join(tab, "PICRUSt2_isoferulic_results.tsv"), sep="\t")
        if iso.robust.any():
            misses.append("PICRUSt2 day17-day25: a hydroxycinnamate feature is robust, but the text says none")

    for f in ("picrust2_nsti", "picrust2_inosine", "picrust2_isoferulic") + tuple(
            f"picrust2_{lv}_{w}" for lv in FILES for w in ("pca", "heatmap", "volcano")):
        if f'src="figures/{f}.png"' not in sec:
            misses.append(f"PICRUSt2 {report}: {f} not shown")
    order = [sec.find(x) for x in ("PICRUSt2 predicts which genes", "Prediction coverage</h3>", "picrust2_nsti.png",
                                   "Whole-profile tests</h3>", "picrust2_MetaCyc_pca.png", "Per-feature tests</h3>",
                                   "<strong>MaAsLin3.</strong>", "MaAsLin3 calls per row", "Robust calls follow",
                                   "picrust2_MetaCyc_heatmap.png", "picrust2_MetaCyc_volcano.png", "Inosine and isoferulic acid</h3>",
                                   "picrust2_inosine.png", "Three limits apply", "<h3>Findings</h3>")]
    if -1 in order or order != sorted(order):
        misses.append(f"PICRUSt2 {report}: section out of order or missing a part: {order}")
    if not (html.index("<h2>Summary</h2>") > p0 > html.index("<h2>Differential abundance</h2>")):
        misses.append(f"PICRUSt2 {report}: section is not between Differential abundance and Summary")
    if sec.count('<div class="finding') != 4:
        misses.append(f"PICRUSt2 {report}: section should close on four findings")
    if re.search(r"LinDA|ggpicrust2 's LinDA", text):
        misses.append(f"PICRUSt2 {report}: section mentions LinDA")
    return misses
