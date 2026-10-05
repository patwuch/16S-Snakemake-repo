#!/usr/bin/env python3
"""Interaction estimates from ANCOM-BC2 against ALDEx3, one point per taxon.

Reads tables/DA_interaction_intersection_Greengenes2.tsv (aldex3_intersection.py).
Both axes are the DietSCFAs:RadiationIR coefficient in log2. Point shade is the
number of mice with a zero count, which is what separates the two methods:
ANCOM-BC2 drops a zero as a missing log, ALDEx3 keeps it. Rings mark q < 0.05
in either method; ringed taxa are named with their zero count.

Usage, inside the qiime2 container:
    python3 plot_aldex3_concordance.py [--genus-names-only]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
OUT_DIR = os.path.join(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
DB = "Greengenes2"
ALPHA = 0.05
MAIN_GAMMA = "0.5"

TEXT_INK = "#2b2b28"
MUTED_INK = "#6b6b66"
# Sequential by zero load: darkest = no zeros, the taxa both methods treat alike.
ZERO_INK = ["#1f3a5f", "#6f8fb3", "#c3cfdd"]
RULE_INK = "#c9c9c4"
ALDEX_HUE = "#377EB8"


def genus_or_above(taxon):
    named = [p.strip() for p in str(taxon).split(";") if p.strip() and not p.strip().endswith("__")]
    return bool(named) and not named[-1].startswith("s__")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genus-names-only", action="store_true",
                    help="only taxa named no further than genus; same fits, same q-values")
    args = ap.parse_args()

    df = pd.read_csv(os.path.join(OUT_DIR, "tables", f"DA_interaction_intersection_{DB}.tsv"), sep="\t")
    if args.genus_names_only:
        df = df[df["taxon"].map(genus_or_above)].reset_index(drop=True)
    x = df["ancombc2_log2fc"]
    y = df[f"aldex3_log2fc_gamma{MAIN_GAMMA}"]
    anc = df["ancombc2_q"] < ALPHA
    ald = df[f"aldex3_q_gamma{MAIN_GAMMA}"] < ALPHA
    both = df["passes_both"]
    zeros = df["n_zero_mice"]
    zclass = {
        "in all 32 mice": (zeros == 0, ZERO_INK[0], 30),
        "zero in 1–8 mice": ((zeros > 0) & (zeros <= 8), ZERO_INK[1], 22),
        "zero in more than 8 mice": (zeros > 8, ZERO_INK[2], 22),
    }

    fig, ax = plt.subplots(figsize=(7.2, 6.8))
    lim = max(abs(x).max(), abs(y).max()) * 1.12
    ax.axhline(0, color=RULE_INK, lw=0.8, zorder=0)
    ax.axvline(0, color=RULE_INK, lw=0.8, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color=RULE_INK, lw=0.8, ls="--", zorder=0)
    for sel, colour, size in zclass.values():
        ax.scatter(x[sel], y[sel], s=size, c=colour, lw=0, zorder=2)
    ax.scatter(x[anc], y[anc], s=95, facecolors="none", edgecolors=TEXT_INK, lw=1.2, zorder=3)
    ax.scatter(x[ald], y[ald], s=95, facecolors="none", edgecolors=ALDEX_HUE, lw=1.2, zorder=3)
    # Labels sit right of their point; nudge apart any two that would overlap.
    named = sorted(df.index[anc | ald], key=lambda i: y[i])
    placed = []
    for i in named:
        ty = y[i]
        for px, py in placed:
            if abs(x[i] - px) < 2.5 and abs(ty - py) < 0.32:
                ty = py + 0.32
        placed.append((x[i], ty))
        ax.annotate(f"{df.at[i, 'label']} ({int(zeros[i])})", (x[i], y[i]), xytext=(x[i] + 0.22, ty),
                    textcoords="data", fontsize=7, color=TEXT_INK, va="center",
                    arrowprops=dict(arrowstyle="-", color=RULE_INK, lw=0.6, shrinkA=0, shrinkB=5))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("ANCOM-BC2 interaction (log2)", color=TEXT_INK, fontsize=9)
    ax.set_ylabel(f"ALDEx3 interaction (log2, CLR scale model, gamma = {MAIN_GAMMA})",
                  color=TEXT_INK, fontsize=9)
    ax.tick_params(colors=MUTED_INK, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color(RULE_INK)
    view = "genus-level names only · " if args.genus_names_only else ""
    ax.set_title(f"Interaction estimates, ANCOM-BC2 against ALDEx3\n"
                 f"{DB} · {view}{len(df)} taxa\n"
                 f"dashed line: equal estimates · (n) = mice with a zero count", fontsize=10, color=TEXT_INK)
    handles = [Line2D([], [], marker="o", ls="none", ms=6 if k.startswith("in") else 5, mfc=c, mec=c,
                      label=f"{k} ({int(sel.sum())})") for k, (sel, c, _) in zclass.items()]
    handles += [
        Line2D([], [], marker="o", ls="none", ms=9, mfc="none", mec=TEXT_INK, mew=1.2,
               label=f"ANCOM-BC2 q < {ALPHA} ({int(anc.sum())})"),
        Line2D([], [], marker="o", ls="none", ms=9, mfc="none", mec=ALDEX_HUE, mew=1.2,
               label=f"ALDEx3 q < {ALPHA} ({int(ald.sum())})"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=7.5, labelcolor=MUTED_INK)
    fig.tight_layout()

    stem = f"DA_interaction_concordance_{DB}" + ("_gplus" if args.genus_names_only else "")
    fig_dir = os.path.join(OUT_DIR, "figures")
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(fig_dir, f"{stem}.{ext}"), dpi=200)
    print(f"Wrote figures/{stem}.png/.svg  ({len(df)} taxa; ANCOM-BC2 {int(anc.sum())}, "
          f"ALDEx3 {int(ald.sum())}, both {int(both.sum())})")


if __name__ == "__main__":
    main()
