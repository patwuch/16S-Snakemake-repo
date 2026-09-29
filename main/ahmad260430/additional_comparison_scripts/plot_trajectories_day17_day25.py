#!/usr/bin/env python3
"""Per-mouse Day17 -> Day25 slopes for the taxa ANCOM-BC2 confirmed, both arms.

Two timepoints give a line segment per mouse, not a curve, so nothing here is
interpolated. What the figure adds over the effect-size barplot is the thing a
single fold change throws away: whether four mice agreed, and how the other arm
behaved on the same taxa.

Why the y axis is not relative abundance. Relative abundance is compositional.
Muribaculaceae giving up forty points mechanically inflates every other taxon in
the same sample, so a rising line can be an artifact of something else falling.
`bias_correct_log_table` is ANCOM-BC2's per-sample log abundance with the
estimated sampling fraction removed, which is the quantity the model actually
compares. The `res` table's `lfc` is bias-corrected too but is one number per
taxon, so it cannot carry a per-mouse view; this can.

The caveat that buys: each arm is a separate fit, so its sampling fraction is
estimated independently. Slopes are comparable across arms, absolute levels only
roughly so.

Usage (from anywhere):
    apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
      /srv/projects/16S-Snakemake-repo/containers/qiime2-2025.10-amplicon-core.sif \
      python3 plot_trajectories_day17_day25.py [--db Greengenes2] [--markers-from IFA_IR]
"""
import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

THIS_DIR = Path(__file__).resolve().parent
STUDY_DIR = THIS_DIR.parent
OUT_DIR = STUDY_DIR / "comparison-day17-day25"

ARMS = {"IFA": "IFA", "IFA_IR": "IFA+IR"}
# The study's own PrimaryColor assignments, as metadata.tsv carries them and the
# factorial report already uses. Colour encodes the ARM here rather than the
# direction of change: the slope and the panel ordering already show direction,
# so colour is free to carry the comparison the figure exists for.
ARM_COLORS = {"IFA": "#90BFF9", "IFA_IR": "#0000FF"}
ARM_DX = {"IFA": -0.055, "IFA_IR": 0.055}   # dodge, so eight lines stay legible
DAYS = ["Day17", "Day25"]
ALPHA = 0.05
NCOL = 4


def mouse_id(sid):
    m = re.match(r"^E\d+_(\d+)\.(\d+)$", sid) or re.match(r"^(\d+)\.(\d+)colon$", sid)
    if not m:
        raise ValueError(f"unparsed sample id: {sid}")
    return f"{m.group(1)}_{m.group(2)}"


def short_label(taxon_path):
    """Deepest rank the classifier actually named.

    Greengenes2 emits a bare `s__` where it resolved no species, so the deepest
    token is often empty. Walking up reports the real call rather than a blank.
    """
    for p in reversed([p.strip() for p in taxon_path.split(";")]):
        name = p.split("__", 1)[-1] if "__" in p else p
        if name:
            return name.replace("_", " "), (p.split("__", 1)[0] if "__" in p else "")
    return taxon_path, ""


def table(db, arm, suffix=""):
    return OUT_DIR / "tables" / f"ANCOMBC2_day17day25_{arm}_{db}{suffix}.tsv"


def load_markers(db, arm):
    f = table(db, arm)
    if not f.is_file():
        raise SystemExit(f"no ANCOM-BC2 table at {f}; run run_ancombc2_day17_day25.R first")
    q, l, ss = "q_TimepointDay25", "lfc_TimepointDay25", "passed_ss_TimepointDay25"
    keep = [{"taxon": r["taxon"], "q": float(r[q]), "log2fc": float(r[l]) / math.log(2)}
            for r in csv.DictReader(open(f), delimiter="\t")
            if r[q] not in ("", "NA") and r[l] not in ("", "NA")
            and float(r[q]) < ALPHA and r.get(ss) in ("TRUE", "True")]
    # Strongest Day25 enrichment first, strongest Day17 enrichment last.
    keep.sort(key=lambda m: -m["log2fc"])
    return keep


def load_bias_corrected(db, arm, meta):
    """taxon -> {(mouse, day): value}. Missing cells are zeros the log drops."""
    f = table(db, arm, "_biascorrected")
    if not f.is_file():
        raise SystemExit(f"no bias-corrected table at {f}; rerun run_ancombc2_day17_day25.R")
    out = {}
    for r in csv.DictReader(open(f), delimiter="\t"):
        vals = {}
        for s, v in r.items():
            if s == "taxon" or s not in meta or v in ("", "NA"):
                continue
            vals[(mouse_id(s), meta[s]["Timepoint"])] = float(v)
        out[r["taxon"]] = vals
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="Greengenes2", choices=["Silva138", "Greengenes2"])
    ap.add_argument("--markers-from", default="IFA_IR", choices=sorted(ARMS))
    args = ap.parse_args()

    figs = OUT_DIR / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    stem = f"ANCOMBC2_day17day25_botharms_{args.db}_trajectories"

    markers = load_markers(args.db, args.markers_from)
    if not markers:
        # Same contract the effect-size barplot honours when nothing is
        # significant: leave the declared outputs behind, empty.
        print(f"No confirmed taxa in {args.markers_from}/{args.db}; writing empty figures.")
        # Truncate, not touch: touch would leave an earlier run's figure in place
        # under a new timestamp, showing taxa this fit no longer confirms.
        for ext in ("svg", "png"):
            (figs / f"{stem}.{ext}").write_bytes(b"")
        return

    meta = {r["sampleid"]: r for r in
            csv.DictReader(open(STUDY_DIR / "metadata.tsv"), delimiter="\t")}
    data = {a: load_bias_corrected(args.db, a, meta) for a in ARMS}
    mice = {a: sorted({mouse_id(s) for s in meta if meta[s]["Group"] == ARMS[a]})
            for a in ARMS}

    vals = [v for a in ARMS for m in markers for v in data[a].get(m["taxon"], {}).values()]
    pad = (max(vals) - min(vals)) * 0.08
    ylim = (min(vals) - pad, max(vals) + pad * 2.4)
    print(f"{len(markers)} taxa from {args.markers_from}/{args.db}; "
          f"shared y {min(vals):.2f} to {max(vals):.2f} (bias-corrected log abundance)")

    nrow = math.ceil(len(markers) / NCOL)
    fig, axes = plt.subplots(nrow, NCOL, figsize=(3.1 * NCOL, 2.6 * nrow),
                             sharex=True, sharey=True)
    axes = axes.ravel()

    for ax, m in zip(axes, markers):
        notes = []
        for arm in ARMS:
            series = data[arm].get(m["taxon"], {})
            agree = total = 0
            for mouse in mice[arm]:
                ys = [series.get((mouse, d)) for d in DAYS]
                xs = [x + ARM_DX[arm] for x in (0, 1)]
                if None not in ys:
                    total += 1
                    if (ys[1] > ys[0]) == (m["log2fc"] > 0):
                        agree += 1
                    ax.plot(xs, ys, color=ARM_COLORS[arm], lw=1.25, alpha=0.8, zorder=2)
                # A missing cell is a zero the log cannot take. Draw the endpoint
                # that exists, hollow, and no segment, rather than inventing one.
                for x, y in zip(xs, ys):
                    if y is None:
                        continue
                    solid = None not in ys
                    ax.plot(x, y, marker="o", ms=4.5, zorder=3,
                            color=ARM_COLORS[arm] if solid else "white",
                            markeredgecolor=ARM_COLORS[arm], markeredgewidth=1.2)
            notes.append((arm, agree, total))

        name, rank = short_label(m["taxon"])
        ax.set_title(name, fontsize=8, pad=3,
                     style="italic" if rank in ("g", "s") else "normal")
        ax.annotate(f"{m['log2fc']:+.2f} log2   q={m['q']:.3f}",
                    xy=(0.5, 0.075), xycoords="axes fraction", ha="center",
                    fontsize=6.5, color="#555555")
        txt = "   ".join(f"{ARMS[a]} {g}/{t}" for a, g, t in notes)
        ax.annotate(txt, xy=(0.5, 0.012), xycoords="axes fraction", ha="center",
                    fontsize=6.5,
                    color="#555555" if all(g == t for _, g, t in notes) else "#9b2226")
        ax.axhline(0, color="#dddddd", lw=0.6, zorder=1)
        ax.set_ylim(*ylim)
        ax.set_xlim(-0.4, 1.4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(DAYS, fontsize=7.5)
        ax.tick_params(axis="y", labelsize=7)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    for ax in axes[len(markers):]:
        ax.set_visible(False)

    fig.suptitle(f"Per-mouse Day17 to Day25, taxa confirmed in {ARMS[args.markers_from]}",
                 fontsize=12.5, y=0.995)
    fig.text(0.5, 0.967,
             f"{args.db} · bias-corrected log abundance, not relative abundance · "
             f"n = 4 mice per arm · counts are mice moving the way the fit says",
             ha="center", fontsize=8.5, color="gray")
    fig.supylabel("Bias-corrected log abundance", fontsize=9)
    fig.legend(handles=[
        Line2D([], [], color=ARM_COLORS["IFA_IR"], lw=1.8, marker="o", ms=5,
               label="IFA+IR (markers from this arm)"),
        Line2D([], [], color=ARM_COLORS["IFA"], lw=1.8, marker="o", ms=5, label="IFA"),
        Line2D([], [], color="#777777", lw=0, marker="o", ms=5, markerfacecolor="white",
               markeredgecolor="#777777", label="other timepoint is zero, no slope drawn"),
    ], loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.004))
    fig.tight_layout(rect=(0.02, 0.028, 1, 0.957))

    for ext, kw in (("svg", {}), ("png", {"dpi": 300})):
        p = figs / f"{stem}.{ext}"
        fig.savefig(p, format=ext, bbox_inches="tight", **kw)
        print(f"  wrote {p}")
    plt.close(fig)


if __name__ == "__main__":
    main()
