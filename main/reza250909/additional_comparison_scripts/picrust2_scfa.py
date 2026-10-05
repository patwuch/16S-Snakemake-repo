#!/usr/bin/env python3
"""Short-chain fatty acid and formate views of the PICRUSt2 LinDA results (reza250909).

For each metabolite (acetate, formate, propionate, butyrate), a curated list of
features in each database (MetaCyc, KO, EC), grouped by role, is read from the
LinDA tables written by run_ggpicrust2_linda.R. The lists were fixed from the
feature annotations alone, before looking at any result.
  * features predicted in >= 16 of 32 mice: the >= 16-mice refit, plus the refits
    without IR mouse C2 (SCFAs+IR vs IR) and without SCFAs mouse N3.2 (SCFAs+IR vs SCFAs);
  * features predicted in fewer mice: the fit with every feature, no one-mouse refits.
    Acetate keeps its own role for these (all nine-mouse features); for the other
    metabolites they are moved to a closing role, their functional role kept in the label.
A feature that serves more than one metabolite is listed under each and says so.

Supersedes picrust2_acetate.py, whose acetate outputs it reproduces unchanged.

Writes, for each metabolite <m>
  tables/PICRUSt2_<m>_features.tsv   the curated list (database, feature, role, label)
  tables/PICRUSt2_<m>_results.tsv    one row per feature x column
  figures/picrust2_<database>_<m>.{png,svg}

Usage, inside the qiime2 container:
    python3 picrust2_scfa.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
TAB = os.path.join(OUT_DIR, "tables")
FIG = os.path.join(OUT_DIR, "figures")

PTA_ACK = "Acetate from acetyl-CoA (Pta–AckA)"
OTHER = "Other acetate-forming enzymes"
WL = "Acetogenesis (Wood–Ljungdahl)"
USE = "Acetate activation, CoA transfer, uptake"
FERM = "Fermentations yielding acetate"
SHARED = "Shares Wood–Ljungdahl enzymes"
SPARSE = "Glyoxylate shunt and acetate import (9 mice)"
FEATURES = [
    ("MetaCyc", "PWY-5100", FERM, "pyruvate fermentation to acetate and lactate II"),
    ("MetaCyc", "P461-PWY", FERM, "hexitol fermentation to lactate, formate, ethanol and acetate"),
    ("MetaCyc", "P163-PWY", FERM, "L-lysine fermentation to acetate and butanoate"),
    ("MetaCyc", "CODH-PWY", WL, "reductive acetyl coenzyme A pathway"),
    ("MetaCyc", "METH-ACETATE-PWY", SHARED, "methanogenesis from acetate (no methanogens predicted)"),
    ("MetaCyc", "PWY-7254", SPARSE, "TCA cycle VII (acetate-producers)"),
    ("MetaCyc", "GLYOXYLATE-BYPASS", SPARSE, "glyoxylate cycle"),
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
    ("KO", "ko:K14393", SPARSE, "actP cation/acetate symporter"),
    ("KO", "ko:K01637", SPARSE, "aceA isocitrate lyase"),
    ("KO", "ko:K01638", SPARSE, "aceB malate synthase"),
    ("KO", "ko:K19709", SPARSE, "ydiF acetate CoA-transferase"),
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
    ("EC", "4.1.3.1", SPARSE, "isocitrate lyase"),
    ("EC", "2.3.3.9", SPARSE, "malate synthase"),
]
SPARSE_LATE = "Predicted in fewer than 16 mice (fit with every feature)"
ALSO_AC = " (also in the acetate panel)"

FORMATE = [
    ("MetaCyc", "FERMENTATION-PWY", "Fermentations yielding formate", "mixed acid fermentation"),
    ("MetaCyc", "P461-PWY", "Fermentations yielding formate", "hexitol fermentation to lactate, formate, ethanol and acetate" + ALSO_AC),
    ("MetaCyc", "1CMET2-PWY", "One-carbon use and release", "N10-formyl-tetrahydrofolate biosynthesis"),
    ("KO", "ko:K00656", "Formate from pyruvate (pyruvate formate-lyase)", "pflB formate C-acetyltransferase"),
    ("KO", "ko:K04069", "Formate from pyruvate (pyruvate formate-lyase)", "pflA pyruvate formate-lyase activating enzyme"),
    ("KO", "ko:K04070", "Formate from pyruvate (pyruvate formate-lyase)", "pflX putative pyruvate formate-lyase activating enzyme"),
    ("KO", "ko:K01577", "Formate from oxalate", "oxc oxalyl-CoA decarboxylase"),
    ("KO", "ko:K07749", "Formate from oxalate", "frc formyl-CoA transferase"),
    ("KO", "ko:K18702", "Formate from oxalate", "uctC CoA:oxalate CoA-transferase"),
    ("KO", "ko:K00123", "Formate oxidation (formate dehydrogenases)", "fdoG/fdhF formate dehydrogenase major subunit"),
    ("KO", "ko:K00124", "Formate oxidation (formate dehydrogenases)", "fdoH formate dehydrogenase iron-sulfur subunit"),
    ("KO", "ko:K00127", "Formate oxidation (formate dehydrogenases)", "fdoI formate dehydrogenase gamma subunit"),
    ("KO", "ko:K05299", "Formate oxidation (formate dehydrogenases)", "fdhA formate dehydrogenase (NADP+) alpha"),
    ("KO", "ko:K15022", "Formate oxidation (formate dehydrogenases)", "fdhB formate dehydrogenase (NADP+) beta"),
    ("KO", "ko:K22339", "Formate oxidation (formate dehydrogenases)", "hylB formate dehydrogenase (NAD+, ferredoxin) B"),
    ("KO", "ko:K22340", "Formate oxidation (formate dehydrogenases)", "hylC formate dehydrogenase (NAD+, ferredoxin) C"),
    ("KO", "ko:K22341", "Formate oxidation (formate dehydrogenases)", "fdhF2 formate dehydrogenase (NAD+, ferredoxin)"),
    ("KO", "ko:K08348", "Formate oxidation (formate dehydrogenases)", "fdnG formate dehydrogenase-N alpha"),
    ("KO", "ko:K08349", "Formate oxidation (formate dehydrogenases)", "fdnH formate dehydrogenase-N beta"),
    ("KO", "ko:K08350", "Formate oxidation (formate dehydrogenases)", "fdnI formate dehydrogenase-N gamma"),
    ("KO", "ko:K22015", "Formate oxidation (formate dehydrogenases)", "fdhF formate dehydrogenase (hydrogenase)"),
    ("KO", "ko:K15827", "Formate hydrogenlyase", "hycB formate hydrogenlyase subunit 2"),
    ("KO", "ko:K15828", "Formate hydrogenlyase", "hycC formate hydrogenlyase subunit 3"),
    ("KO", "ko:K15829", "Formate hydrogenlyase", "hycD formate hydrogenlyase subunit 4"),
    ("KO", "ko:K15830", "Formate hydrogenlyase", "hycE formate hydrogenlyase subunit 5"),
    ("KO", "ko:K15831", "Formate hydrogenlyase", "hycF formate hydrogenlyase subunit 6"),
    ("KO", "ko:K15832", "Formate hydrogenlyase", "hycG formate hydrogenlyase subunit 7"),
    ("KO", "ko:K01938", "One-carbon use and release", "fhs formate–tetrahydrofolate ligase" + ALSO_AC),
    ("KO", "ko:K00288", "One-carbon use and release", "MTHFD trifunctional C1-tetrahydrofolate synthase"),
    ("KO", "ko:K01433", "One-carbon use and release", "purU formyltetrahydrofolate deformylase"),
    ("KO", "ko:K21993", "Formate transport", "fdhC formate transporter"),
    ("KO", "ko:K06212", "Formate transport", "focA formate channel"),
    ("KO", "ko:K08177", "Formate transport", "oxlT oxalate/formate antiporter"),
    ("KO", "ko:K21990", "Formate transport", "yfdC formate–nitrite transporter family protein"),
    ("EC", "2.3.1.54", "Formate from pyruvate (pyruvate formate-lyase)", "formate C-acetyltransferase"),
    ("EC", "1.97.1.4", "Formate from pyruvate (pyruvate formate-lyase)", "pyruvate formate-lyase activating enzyme"),
    ("EC", "4.1.1.8", "Formate from oxalate", "oxalyl-CoA decarboxylase"),
    ("EC", "2.8.3.16", "Formate from oxalate", "formyl-CoA transferase"),
    ("EC", "2.8.3.19", "Formate from oxalate", "CoA:oxalate CoA-transferase"),
    ("EC", "1.17.1.9", "Formate oxidation (formate dehydrogenases)", "formate dehydrogenase"),
    ("EC", "1.17.1.10", "Formate oxidation (formate dehydrogenases)", "formate dehydrogenase (NADP+)"),
    ("EC", "1.17.1.11", "Formate oxidation (formate dehydrogenases)", "formate dehydrogenase (NAD+, ferredoxin)"),
    ("EC", "1.17.5.3", "Formate oxidation (formate dehydrogenases)", "formate dehydrogenase-N"),
    ("EC", "6.3.4.3", "One-carbon use and release", "formate–tetrahydrofolate ligase" + ALSO_AC),
    ("EC", "3.5.1.10", "One-carbon use and release", "formyltetrahydrofolate deformylase"),
]

PROPIONATE = [
    ("MetaCyc", "P108-PWY", "Succinate pathway", "pyruvate fermentation to propanoate I"),
    ("MetaCyc", "PWY-5494", "Acrylate pathway", "pyruvate fermentation to propanoate II (acrylate pathway)"),
    ("MetaCyc", "PWY-7013", "Propanediol pathway", "L-1,2-propanediol degradation"),
    ("MetaCyc", "PWY-5130", "From amino acids", "2-oxobutanoate degradation I"),
    ("MetaCyc", "PWY0-42", "Propionate use (methylcitrate cycle)", "2-methylcitrate cycle I"),
    ("MetaCyc", "PWY-5747", "Propionate use (methylcitrate cycle)", "2-methylcitrate cycle II"),
    ("KO", "ko:K01847", "Succinate pathway", "MUT methylmalonyl-CoA mutase"),
    ("KO", "ko:K01848", "Succinate pathway", "mcmA1 methylmalonyl-CoA mutase, N-terminal domain"),
    ("KO", "ko:K01849", "Succinate pathway", "mcmA2 methylmalonyl-CoA mutase, C-terminal domain"),
    ("KO", "ko:K05606", "Succinate pathway", "epi methylmalonyl-CoA epimerase"),
    ("KO", "ko:K11264", "Succinate pathway", "scpB/mmcD methylmalonyl-CoA decarboxylase"),
    ("KO", "ko:K03416", "Succinate pathway", "methylmalonyl-CoA carboxyltransferase 5S subunit"),
    ("KO", "ko:K20626", "Acrylate pathway", "lcdA lactoyl-CoA dehydratase alpha"),
    ("KO", "ko:K20627", "Acrylate pathway", "lcdB lactoyl-CoA dehydratase beta"),
    ("KO", "ko:K01026", "Acrylate pathway", "pct propionate CoA-transferase"),
    ("KO", "ko:K01699", "Propanediol pathway", "pduC propanediol dehydratase large subunit"),
    ("KO", "ko:K13919", "Propanediol pathway", "pduD propanediol dehydratase medium subunit"),
    ("KO", "ko:K13920", "Propanediol pathway", "pduE propanediol dehydratase small subunit"),
    ("KO", "ko:K13922", "Propanediol pathway", "pduP propionaldehyde dehydrogenase"),
    ("KO", "ko:K13921", "Propanediol pathway", "pduQ 1-propanol dehydrogenase"),
    ("KO", "ko:K13923", "Propanediol pathway", "pduL phosphate propanoyltransferase"),
    ("KO", "ko:K00932", "Propanediol pathway", "tdcD propionate kinase"),
    ("KO", "ko:K01965", "Propionate use (methylcitrate cycle)", "pccA propionyl-CoA carboxylase alpha"),
    ("EC", "5.4.99.2", "Succinate pathway", "methylmalonyl-CoA mutase"),
    ("EC", "5.1.99.1", "Succinate pathway", "methylmalonyl-CoA epimerase"),
    ("EC", "4.1.1.41", "Succinate pathway", "methylmalonyl-CoA decarboxylase"),
    ("EC", "2.1.3.1", "Succinate pathway", "methylmalonyl-CoA carboxytransferase"),
    ("EC", "4.2.1.54", "Acrylate pathway", "lactoyl-CoA dehydratase"),
    ("EC", "2.8.3.1", "Acrylate pathway", "propionate CoA-transferase"),
    ("EC", "4.2.1.28", "Propanediol pathway", "propanediol dehydratase"),
    ("EC", "1.2.1.87", "Propanediol pathway", "propanal dehydrogenase (CoA-propanoylating)"),
    ("EC", "2.3.1.222", "Propanediol pathway", "phosphate propanoyltransferase"),
    ("EC", "2.7.2.15", "Propanediol pathway", "propionate kinase"),
    ("EC", "6.4.1.3", "Propionate use (methylcitrate cycle)", "propionyl-CoA carboxylase"),
]

BUTYRATE = [
    ("MetaCyc", "P163-PWY", "Lysine pathway", "L-lysine fermentation to acetate and butanoate" + ALSO_AC),
    ("MetaCyc", "PWY-5677", "4-aminobutyrate and succinate pathways", "succinate fermentation to butanoate"),
    ("MetaCyc", "PWY-5022", "4-aminobutyrate and succinate pathways", "4-aminobutanoate degradation V"),
    ("KO", "ko:K00626", "Acetyl-CoA to butyryl-CoA", "thl/atoB acetyl-CoA C-acetyltransferase"),
    ("KO", "ko:K00074", "Acetyl-CoA to butyryl-CoA", "hbd 3-hydroxybutyryl-CoA dehydrogenase"),
    ("KO", "ko:K01715", "Acetyl-CoA to butyryl-CoA", "crt enoyl-CoA hydratase"),
    ("KO", "ko:K00248", "Acetyl-CoA to butyryl-CoA", "bcd butyryl-CoA dehydrogenase"),
    ("KO", "ko:K03522", "Acetyl-CoA to butyryl-CoA", "etfA electron transfer flavoprotein alpha"),
    ("KO", "ko:K03521", "Acetyl-CoA to butyryl-CoA", "etfB electron transfer flavoprotein beta"),
    ("KO", "ko:K00634", "Butyrate from butyryl-CoA", "ptb phosphate butyryltransferase"),
    ("KO", "ko:K00929", "Butyrate from butyryl-CoA", "buk butyrate kinase"),
    ("KO", "ko:K01034", "Butyrate from butyryl-CoA", "atoD acetate/acetoacetate CoA-transferase alpha" + ALSO_AC),
    ("KO", "ko:K01035", "Butyrate from butyryl-CoA", "atoA acetate/acetoacetate CoA-transferase beta" + ALSO_AC),
    ("KO", "ko:K19709", "Butyrate from butyryl-CoA", "ydiF acetate CoA-transferase" + ALSO_AC),
    ("KO", "ko:K01843", "Lysine pathway", "kamA lysine 2,3-aminomutase"),
    ("KO", "ko:K01844", "Lysine pathway", "kamD beta-lysine 5,6-aminomutase alpha"),
    ("KO", "ko:K18011", "Lysine pathway", "kamE beta-lysine 5,6-aminomutase beta"),
    ("KO", "ko:K18012", "Lysine pathway", "kdd 3,5-diaminohexanoate dehydrogenase"),
    ("KO", "ko:K18013", "Lysine pathway", "kce 3-keto-5-aminohexanoate cleavage enzyme"),
    ("KO", "ko:K18014", "Lysine pathway", "kal 3-aminobutyryl-CoA ammonia-lyase"),
    ("KO", "ko:K01039", "Glutarate pathway", "gctA glutaconate CoA-transferase A"),
    ("KO", "ko:K01040", "Glutarate pathway", "gctB glutaconate CoA-transferase B"),
    ("KO", "ko:K01615", "Glutarate pathway", "gcdA glutaconyl-CoA decarboxylase alpha"),
    ("KO", "ko:K14534", "4-aminobutyrate and succinate pathways", "abfD 4-hydroxybutyryl-CoA dehydratase"),
    ("KO", "ko:K18122", "4-aminobutyrate and succinate pathways", "abfT 4-hydroxybutyrate CoA-transferase"),
    ("EC", "2.3.1.9", "Acetyl-CoA to butyryl-CoA", "acetyl-CoA C-acetyltransferase"),
    ("EC", "1.1.1.157", "Acetyl-CoA to butyryl-CoA", "3-hydroxybutyryl-CoA dehydrogenase"),
    ("EC", "4.2.1.55", "Acetyl-CoA to butyryl-CoA", "3-hydroxybutyryl-CoA dehydratase"),
    ("EC", "4.2.1.17", "Acetyl-CoA to butyryl-CoA", "enoyl-CoA hydratase"),
    ("EC", "1.3.8.1", "Acetyl-CoA to butyryl-CoA", "short-chain acyl-CoA dehydrogenase"),
    ("EC", "1.3.1.44", "Acetyl-CoA to butyryl-CoA", "trans-2-enoyl-CoA reductase (NAD+)"),
    ("EC", "2.3.1.19", "Butyrate from butyryl-CoA", "phosphate butyryltransferase"),
    ("EC", "2.7.2.7", "Butyrate from butyryl-CoA", "butyrate kinase"),
    ("EC", "2.8.3.8", "Butyrate from butyryl-CoA", "acetate CoA-transferase (butyryl-CoA:acetate)" + ALSO_AC),
    ("EC", "5.4.3.2", "Lysine pathway", "lysine 2,3-aminomutase"),
    ("EC", "2.8.3.12", "Glutarate pathway", "glutaconate CoA-transferase"),
    ("EC", "4.2.1.120", "4-aminobutyrate and succinate pathways", "4-hydroxybutanoyl-CoA dehydratase"),
]

# name -> (features, heading word, role kept for sparse features or None to move them to SPARSE_LATE)
SETS = {
    "acetate": (FEATURES, "Acetate", SPARSE),
    "formate": (FORMATE, "Formate", None),
    "propionate": (PROPIONATE, "Propionate", None),
    "butyrate": (BUTYRATE, "Butyrate", None),
}
COLUMNS = [("A", "SCFAs vs Control", ""), ("B", "SCFAs+IR vs IR", ""), ("C", "IR vs Control", ""),
           ("D", "SCFAs+IR vs SCFAs", ""), ("B′", "SCFAs+IR vs IR", "_exclC2"), ("D′", "SCFAs+IR vs SCFAs", "_exclN3.2")]
COL_LABEL = {"A": "SCFAs vs\nControl", "B": "SCFAs+IR\nvs IR", "C": "IR vs\nControl", "D": "SCFAs+IR\nvs SCFAs",
             "B′": "SCFAs+IR\nvs IR\nwithout C2", "D′": "SCFAs+IR\nvs SCFAs\nwithout N3.2"}
INK = "#14181d"


def load(db, suffix):
    return pd.read_csv(os.path.join(TAB, f"PICRUSt2_{db}_linda{suffix}.tsv"), sep="\t", quoting=3)


def run(name, features, heading, sparse_role):
    feats = pd.DataFrame(features, columns=["database", "feature", "role", "label"])
    if feats.duplicated(["database", "feature"]).any():
        raise SystemExit(f"{name}: duplicate features")
    if sparse_role is None:
        # Features predicted in fewer than 16 mice come from a different fit, so
        # they close the list under their own heading, their role kept in the label.
        for db in ("MetaCyc", "KO", "EC"):
            npred = load(db, "").drop_duplicates("feature").set_index("feature").n_predicted
            sel = feats.database == db
            missing = set(feats.feature[sel]) - set(npred.index)
            if missing:
                raise SystemExit(f"{name} {db}: not predicted by PICRUSt2: {sorted(missing)}")
            few = sel & feats.feature.map(lambda f: npred.get(f, 99) < 16)
            feats.loc[few, "label"] = [f"{l} [{r.split(' (')[0].lower()}; {int(npred[f])} mice]"
                                       for l, r, f in zip(feats.label[few], feats.role[few], feats.feature[few])]
            feats.loc[few, "role"] = SPARSE_LATE
        feats = pd.concat([feats[feats.role != SPARSE_LATE], feats[feats.role == SPARSE_LATE]], ignore_index=True)
        sparse_role = SPARSE_LATE
    feats.to_csv(os.path.join(TAB, f"PICRUSt2_{name}_features.tsv"), sep="\t", index=False)
    rows = []
    for db in ("MetaCyc", "KO", "EC"):
        allf = load(db, "")
        npred = allf.drop_duplicates("feature").set_index("feature").n_predicted
        fits = {s: load(db, "_minpres16" + s) for s in ("", "_exclC2", "_exclN3.2")}
        for r in feats[feats.database == db].itertuples():
            sparse = npred[r.feature] < 16
            if sparse != (r.role == sparse_role):
                raise SystemExit(f"{r.feature}: predicted in {npred[r.feature]} mice, role {r.role}")
            for col, contrast, suf in COLUMNS:
                src = allf if sparse else fits[suf]
                if sparse and suf:
                    d = None
                else:
                    x = src[(src.feature == r.feature) & (src.contrast == contrast)]
                    d = x.iloc[0] if len(x) else None
                rows.append(dict(database=db, feature=r.feature, role=r.role, label=r.label,
                                 n_predicted=int(npred[r.feature]), column=col, contrast=contrast,
                                 fit=("all features" if sparse else "predicted in >= 16 mice") + (" " + suf[1:] if suf else ""),
                                 log2fc=d.log2fc if d is not None else np.nan, p=d.p if d is not None else np.nan,
                                 q=d.q if d is not None else np.nan))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(TAB, f"PICRUSt2_{name}_results.tsv"), sep="\t", index=False)

    norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
    cmap = plt.get_cmap("RdBu_r")
    for db in ("MetaCyc", "KO", "EC"):
        f_ = feats[feats.database == db]
        # One heading row per role, then its features.
        layout = []
        for ro in dict.fromkeys(f_.role):
            layout.append(("role", ro, None))
            layout += [("feat", r.label, r.feature) for r in f_[f_.role == ro].itertuples()]
        n, m = len(layout), len(COLUMNS)
        fig, ax = plt.subplots(figsize=(11, 0.4 * n + 2.2))
        outlines = []
        for i, (kind, text, fid) in enumerate(layout):
            if kind == "role":
                ax.text(-0.1, i + 0.62, text, ha="right", va="center", fontsize=9, fontweight="bold", color=INK)
                ax.axhline(i + 1, color="#9ca3af", lw=0.6, xmin=0, xmax=1)
                continue
            ax.text(-0.1, i + 0.5, text, ha="right", va="center", fontsize=8.5, color=INK)
            for j, (col, _c, _s) in enumerate(COLUMNS):
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
        # Significance outlines on top, inset so neighbouring cells do not cover them.
        for j, i in outlines:
            ax.add_patch(plt.Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, facecolor="none", edgecolor=INK, lw=1.6, zorder=5))
        ax.axvline(4, color=INK, lw=1.2)
        ax.set_xlim(0, m); ax.set_ylim(n, 0)
        ax.set_xticks(np.arange(m) + 0.5); ax.set_xticklabels([COL_LABEL[c] for c, _, _ in COLUMNS], fontsize=8.5)
        ax.xaxis.tick_top()
        ax.set_yticks([])
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
        cb.set_label("log2 fold change, first-named group against the second (clipped at ±2)", fontsize=8)
        ax.set_title(f"{heading}-related predicted features, {db}: log2 fold change by contrast; outlined and bold where raw p < 0.05",
                     fontsize=9.5, pad=34, loc="left")
        fig.tight_layout()
        for ext, kw in (("svg", {}), ("png", {"dpi": 300})):
            fig.savefig(os.path.join(FIG, f"picrust2_{db}_{name}.{ext}"), bbox_inches="tight", **kw)
        plt.close(fig)
        print(f"wrote picrust2_{db}_{name} ({sum(1 for k, _, _ in layout if k == 'feat')} features)")


def main():
    for name, (features, heading, sparse_role) in SETS.items():
        run(name, features, heading, sparse_role)


if __name__ == "__main__":
    main()
