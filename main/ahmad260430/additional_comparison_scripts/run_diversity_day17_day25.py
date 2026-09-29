#!/usr/bin/env python3
"""Day17-vs-Day25 diversity statistics, held within each treatment arm.

Nothing in the comparison-* folders carries a p-value: plot_diversity.py renders
figures only. This computes the tests those figures imply, plus the paired
versions the design actually supports (same four mice sampled twice).

q-values are Benjamini-Hochberg within each arm, across the metrics the report
reads: the four database-independent ones plus the Greengenes2 tree metrics.
Silva138 rows are still written but carry no q, the same family rule the
factorial analysis uses (its log, "Diversity BH families exclude Silva138").
"""
import itertools, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from qiime2 import Artifact
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, wilcoxon
from skbio.stats.distance import permanova, permdisp
from skbio.stats.distance._permdisp import geomedian_axis_one
from skbio.stats.ordination import pcoa
from statsmodels.stats.multitest import multipletests

STUDY = Path(__file__).resolve().parent.parent
CM = STUDY / "qiime2_artifacts" / "core-metrics-results"
ARMS = {"IFA": "IFA", "IFA+IR": "IFA_IR"}   # display name -> filename slug

ALPHA = {  # metric -> (artifact, is_db_specific)
    "shannon": ("shannon-vector.qza", False),
    "chao1": ("chao1-vector.qza", False),
    "evenness": ("evenness-vector.qza", False),
    "simpson": ("simpson-vector.qza", False),
    "faith_pd": ("{db}-faith-pd-vector.qza", True),
}
BETA = {
    "jaccard": ("jaccard-distance-matrix.qza", False),
    "braycurtis": ("bray-curtis-distance-matrix.qza", False),
    "unweighted_unifrac": ("{db}-unweighted-unifrac-distance-matrix.qza", True),
    "weighted_unifrac": ("{db}-weighted-unifrac-distance-matrix.qza", True),
}
DBS = ["Silva138", "Greengenes2"]
SEED = 42


def mouse_id(sid):
    """E1_1.00 and 1.00colon are the same mouse at two timepoints."""
    m = re.match(r"^E\d+_(\d+)\.(\d+)$", sid) or re.match(r"^(\d+)\.(\d+)colon$", sid)
    if not m:
        raise ValueError(f"unparsed sample id: {sid}")
    return f"{m.group(1)}_{m.group(2)}"


def load(rel, db):
    return Artifact.load(str(CM / rel.format(db=db))).view(
        pd.Series if "vector" in rel else __import__("skbio").DistanceMatrix)


def distances_to_median(dm, grouping):
    """Per-sample distance to its group's spatial median, the quantity permdisp tests.

    Rebuilt from skbio's own steps (eigh PCoA over every axis, geometric median
    per group, Euclidean distance) so the dispersion figure plots exactly what
    permdisp's F is computed from.
    """
    coords = pcoa(dm, method="eigh", number_of_dimensions=0).samples
    out = {}
    for label in sorted(set(grouping)):
        ids = [i for i, g in zip(dm.ids, grouping) if g == label]
        X = coords.loc[ids].values.astype("float64")
        med = np.asarray(geomedian_axis_one(X.T))
        out.update(zip(ids, cdist(X, [med]).ravel()))
    return out


def bh_within_arm(df, pcol, qcol):
    df[qcol] = np.nan
    for arm, idx in df[df.db != "Silva138"].groupby("arm").groups.items():
        ps = df.loc[idx, pcol]
        ok = ps.notna()
        df.loc[ps.index[ok], qcol] = multipletests(ps[ok], method="fdr_bh")[1]


rows_a, rows_b, rows_d = [], [], []
for arm, folder in ARMS.items():
    md = pd.read_csv(STUDY / "comparison-day17-day25" / "tables" /
                     f"metadata_subset_{folder}.tsv", sep="\t", index_col=0)
    md["MouseID"] = [mouse_id(s) for s in md.index]
    d17 = md.index[md.Timepoint == "Day17"].tolist()
    d25 = md.index[md.Timepoint == "Day25"].tolist()
    pairs = sorted(md.MouseID.unique())

    for metric, (rel, db_specific) in ALPHA.items():
        for db in (DBS if db_specific else [None]):
            s = load(rel, db)
            a = np.array([s[x] for x in d17], float)
            b = np.array([s[x] for x in d25], float)
            # paired: order both arrays by mouse
            pa = np.array([s[md.index[(md.MouseID == m) & (md.Timepoint == "Day17")][0]] for m in pairs], float)
            pb = np.array([s[md.index[(md.MouseID == m) & (md.Timepoint == "Day25")][0]] for m in pairs], float)
            u_p = mannwhitneyu(a, b, alternative="two-sided").pvalue
            w_p = wilcoxon(pa, pb).pvalue if not np.allclose(pa, pb) else np.nan
            rows_a.append(dict(arm=arm, metric=metric, db=db or "db-independent",
                               day17_mean=a.mean(), day25_mean=b.mean(),
                               delta=b.mean() - a.mean(),
                               mwu_p=u_p, wilcoxon_paired_p=w_p,
                               n_pairs_same_direction=int(max((pb > pa).sum(), (pb < pa).sum()))))

    for metric, (rel, db_specific) in BETA.items():
        for db in (DBS if db_specific else [None]):
            dm = load(rel, db).filter(md.index.tolist())
            g = md.loc[list(dm.ids), "Timepoint"]
            # skbio 0.6 takes no seed argument; it permutes with numpy's global state.
            np.random.seed(SEED)
            pn = permanova(dm, g, permutations=9999)
            np.random.seed(SEED)
            pd_ = permdisp(dm, g, permutations=9999)
            n_s, n_g = len(dm.ids), g.nunique()
            f = pn["test statistic"]
            r2 = f * (n_g - 1) / (f * (n_g - 1) + (n_s - n_g))
            for sid, dist in distances_to_median(dm, list(g)).items():
                rows_d.append(dict(arm=arm, metric=metric, db=db or "db-independent",
                                   sample=sid, MouseID=md.loc[sid, "MouseID"],
                                   Timepoint=md.loc[sid, "Timepoint"], distance=dist))
            within17 = [dm[i, j] for i, j in itertools.combinations(d17, 2)]
            within25 = [dm[i, j] for i, j in itertools.combinations(d25, 2)]
            between = [dm[i, j] for i in d17 for j in d25]
            rows_b.append(dict(arm=arm, metric=metric, db=db or "db-independent",
                               pseudo_F=f, R2=r2, permanova_p=pn["p-value"],
                               permdisp_p=pd_["p-value"],
                               mean_within_d17=np.mean(within17),
                               mean_within_d25=np.mean(within25),
                               mean_between=np.mean(between)))

A = pd.DataFrame(rows_a)
B = pd.DataFrame(rows_b)
D = pd.DataFrame(rows_d)
bh_within_arm(A, "mwu_p", "mwu_q")
bh_within_arm(A, "wilcoxon_paired_p", "wilcoxon_paired_q")
bh_within_arm(B, "permanova_p", "permanova_q")
bh_within_arm(B, "permdisp_p", "permdisp_q")
out = Path(sys.argv[1])
(out / "tables").mkdir(parents=True, exist_ok=True)
A.to_csv(out / "tables" / "diversity_day17day25_alpha_stats.tsv", sep="\t", index=False, float_format="%.6g")
B.to_csv(out / "tables" / "diversity_day17day25_beta_stats.tsv", sep="\t", index=False, float_format="%.6g")
D.to_csv(out / "tables" / "diversity_day17day25_dispersion_distances.tsv", sep="\t", index=False, float_format="%.6g")
pd.set_option("display.width", 200, "display.max_columns", 50)
print("=== ALPHA ===\n", A.to_string(index=False), sep="")
print("\n=== BETA ===\n", B.to_string(index=False), sep="")
