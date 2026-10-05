#!/usr/bin/env python3
"""Acetate-focused view of the PICRUSt2 LinDA results, for either ahmad260430 report.

The curated list of acetate features in each database (MetaCyc, KO, EC), grouped by
role, is the one used for reza250909, chosen by role before looking at any result.
Values come from the LinDA tables written by run_picrust2_linda.R:
  * features predicted in at least half the report's samples: the refit among those
    features, plus any one-sample refits listed in REFITS;
  * features predicted in fewer samples: the fit with every feature, no refits.
A feature's role is kept; whether it is sparse is read from the data.

Writes, in the report's folder,
  tables/PICRUSt2_acetate_features.tsv   the curated list (database, feature, role, label)
  tables/PICRUSt2_acetate_results.tsv    one row per feature x column
  figures/picrust2_<database>_acetate.{png,svg}

Usage, inside the qiime2 container:
    python3 picrust2_acetate.py factorial|day17-day25
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from picrust2_report import REPORTS

PTA_ACK = "Acetate from acetyl-CoA (Pta–AckA)"
OTHER = "Other acetate-forming enzymes"
WL = "Acetogenesis (Wood–Ljungdahl)"
USE = "Acetate activation, CoA transfer, uptake"
FERM = "Fermentations yielding acetate"
SHARED = "Shares Wood–Ljungdahl enzymes"
GLYOX = "Glyoxylate shunt and acetate import"
SPARSE = "Predicted in fewer than half the samples"
FEATURES = [
    ("MetaCyc", "PWY-5100", FERM, "pyruvate fermentation to acetate and lactate II"),
    ("MetaCyc", "P461-PWY", FERM, "hexitol fermentation to lactate, formate, ethanol and acetate"),
    ("MetaCyc", "P163-PWY", FERM, "L-lysine fermentation to acetate and butanoate"),
    ("MetaCyc", "CODH-PWY", WL, "reductive acetyl coenzyme A pathway"),
    ("MetaCyc", "METH-ACETATE-PWY", SHARED, "methanogenesis from acetate (no methanogens predicted)"),
    ("MetaCyc", "PWY-7254", GLYOX, "TCA cycle VII (acetate-producers)"),
    ("MetaCyc", "GLYOXYLATE-BYPASS", GLYOX, "glyoxylate cycle"),
    ("KO", "ko:K00925", PTA_ACK, "ackA acetate kinase"),
    ("KO", "ko:K00625", PTA_ACK, "pta phosphate acetyltransferase"),
    ("KO", "ko:K13788", PTA_ACK, "pta phosphate acetyltransferase (K13788)"),
    ("KO", "ko:K00156", OTHER, "poxB pyruvate dehydrogenase (quinone)"),
    ("KO", "ko:K00158", OTHER, "spxB pyruvate oxidase"),
    ("KO", "ko:K01067", OTHER, "acetyl-CoA hydrolase"),
    ("KO", "ko:K01905", OTHER, "acdA acetate–CoA ligase (ADP-forming)"),
    ("KO", "ko:K18118", OTHER, "aarC succinyl-CoA:acetate CoA-transferase"),
    ("KO", "ko:K14138", WL, "acsB acetyl-CoA synthase"),
    ("KO", "ko:K00198", WL, "cooS/acsA CO dehydrogenase"),
    ("KO", "ko:K01938", WL, "fhs formate–tetrahydrofolate ligase"),
    ("KO", "ko:K01895", USE, "acs acetyl-CoA synthetase"),
    ("KO", "ko:K01034", USE, "atoD acetate CoA-transferase (alpha)"),
    ("KO", "ko:K01035", USE, "atoA acetate CoA-transferase (beta)"),
    ("KO", "ko:K07034", USE, "acetate uptake transporter"),
    ("KO", "ko:K14393", GLYOX, "actP cation/acetate symporter"),
    ("KO", "ko:K01637", GLYOX, "aceA isocitrate lyase"),
    ("KO", "ko:K01638", GLYOX, "aceB malate synthase"),
    ("KO", "ko:K19709", GLYOX, "ydiF acetate CoA-transferase"),
    ("EC", "2.7.2.1", PTA_ACK, "acetate kinase"),
    ("EC", "2.3.1.8", PTA_ACK, "phosphate acetyltransferase"),
    ("EC", "1.2.5.1", OTHER, "pyruvate dehydrogenase (quinone)"),
    ("EC", "1.2.3.3", OTHER, "pyruvate oxidase"),
    ("EC", "3.1.2.1", OTHER, "acetyl-CoA hydrolase"),
    ("EC", "6.2.1.13", OTHER, "acetate–CoA ligase (ADP-forming)"),
    ("EC", "2.8.3.18", OTHER, "succinyl-CoA:acetate CoA-transferase"),
    ("EC", "2.3.1.169", WL, "CO-methylating acetyl-CoA synthase"),
    ("EC", "1.2.7.4", WL, "anaerobic CO dehydrogenase"),
    ("EC", "6.3.4.3", WL, "formate–tetrahydrofolate ligase"),
    ("EC", "6.2.1.1", USE, "acetate–CoA ligase (acs)"),
    ("EC", "2.8.3.8", USE, "acetate CoA-transferase"),
    ("EC", "4.1.3.1", GLYOX, "isocitrate lyase"),
    ("EC", "2.3.3.9", GLYOX, "malate synthase"),
]
# One-sample refits shown as extra columns: report -> [(column, contrast, exclusion suffix, heading)].
REFITS = {"factorial": [], "day17-day25": []}
INK = "#14181d"


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in REPORTS:
        sys.exit(f"usage: picrust2_acetate.py {'|'.join(REPORTS)}")
    rep = REPORTS[sys.argv[1]]
    tab, fig_dir = os.path.join(rep.out, "tables"), os.path.join(rep.out, "figures")
    load = lambda db, suffix: pd.read_csv(os.path.join(tab, f"PICRUSt2_{db}_linda{suffix}.tsv"), sep="\t", quoting=3)
    columns = [(chr(65 + i), c.name, "", c.name) for i, c in enumerate(rep.contrasts)] + REFITS[rep.key]
    feats = pd.DataFrame(FEATURES, columns=["database", "feature", "role", "label"])
    rows = []
    for db in ("MetaCyc", "KO", "EC"):
        allf = load(db, "")
        npred = allf.drop_duplicates("feature").set_index("feature").n_predicted
        fits = {s: load(db, f"_minpres{rep.half}" + s) for s in {c[2] for c in columns}}
        for r in feats[feats.database == db].itertuples():
            if r.feature not in npred.index:
                continue
            sparse = npred[r.feature] < rep.half
            for col, contrast, suf, _ in columns:
                d = None
                if not (sparse and suf):
                    src = allf if sparse else fits[suf]
                    x = src[(src.feature == r.feature) & (src.contrast == contrast)]
                    d = x.iloc[0] if len(x) else None
                rows.append(dict(database=db, feature=r.feature, role=SPARSE if sparse else r.role, label=r.label,
                                 n_predicted=int(npred[r.feature]), column=col, contrast=contrast,
                                 fit=("all features" if sparse else f"predicted in >= {rep.half} samples") + (" " + suf[1:] if suf else ""),
                                 log2fc=d.log2fc if d is not None else np.nan, p=d.p if d is not None else np.nan,
                                 q=d.q if d is not None else np.nan))
    res = pd.DataFrame(rows)
    feats = feats.merge(res.drop_duplicates(["database", "feature"])[["database", "feature", "role", "n_predicted"]],
                        on=["database", "feature"], suffixes=("_curated", ""))
    feats.to_csv(os.path.join(tab, "PICRUSt2_acetate_features.tsv"), sep="\t", index=False)
    res.to_csv(os.path.join(tab, "PICRUSt2_acetate_results.tsv"), sep="\t", index=False)

    norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
    cmap = plt.get_cmap("RdBu_r")
    for db in ("MetaCyc", "KO", "EC"):
        f_ = res[res.database == db].drop_duplicates("feature")
        layout = []
        for ro in dict.fromkeys(f_.role):
            layout.append(("role", ro, None))
            layout += [("feat", r.label, r.feature) for r in f_[f_.role == ro].itertuples()]
        n, m = len(layout), len(columns)
        fig, ax = plt.subplots(figsize=(3.2 + 1.25 * m, 0.4 * n + 2.4))
        outlines = []
        for i, (kind, text, fid) in enumerate(layout):
            if kind == "role":
                ax.text(-0.1, i + 0.62, text, ha="right", va="center", fontsize=9, fontweight="bold", color=INK)
                ax.axhline(i + 1, color="#9ca3af", lw=0.6, xmin=0, xmax=1)
                continue
            ax.text(-0.1, i + 0.5, text, ha="right", va="center", fontsize=8.5, color=INK)
            for j, (col, _c, _s, _h) in enumerate(columns):
                d = res[(res.database == db) & (res.feature == fid) & (res.column == col)].iloc[0]
                if np.isnan(d.log2fc):
                    ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#f3f4f6", edgecolor="white", lw=2))
                    ax.text(j + 0.5, i + 0.5, "—", ha="center", va="center", fontsize=8, color="#9ca3af")
                    continue
                sig = d.p < 0.05
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=cmap(norm(np.clip(d.log2fc, -2, 2))),
                                           edgecolor="white", lw=2))
                ax.text(j + 0.5, i + 0.5, f"{d.log2fc:+.2f}", ha="center", va="center", fontsize=8,
                        color="white" if abs(d.log2fc) > 1.3 else INK, fontweight="bold" if sig else "normal")
                if sig:
                    outlines.append((j, i))
        for j, i in outlines:
            ax.add_patch(plt.Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, facecolor="none", edgecolor=INK, lw=1.6, zorder=5))
        if REFITS[rep.key]:
            ax.axvline(len(rep.contrasts), color=INK, lw=1.2)
        ax.set_xlim(0, m); ax.set_ylim(n, 0)
        ax.set_xticks(np.arange(m) + 0.5)
        ax.set_xticklabels([f"{col}\n" + h.replace(" vs ", "\nvs ") for col, _c, _s, h in columns], fontsize=8)
        ax.xaxis.tick_top()
        ax.set_yticks([])
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
        cb.set_label("log2 fold change, first-named against the second (clipped at ±2)", fontsize=8)
        ax.set_title(f"Acetate-related predicted features, {db}: log2 fold change by contrast; outlined and bold where raw p < 0.05",
                     fontsize=9.5, pad=46, loc="left")
        fig.tight_layout()
        for ext, kw in (("svg", {}), ("png", {"dpi": 300})):
            fig.savefig(os.path.join(fig_dir, f"picrust2_{db}_acetate.{ext}"), bbox_inches="tight", **kw)
        plt.close(fig)
        print(f"wrote picrust2_{db}_acetate ({sum(1 for k, _, _ in layout if k == 'feat')} features)")


if __name__ == "__main__":
    main()
