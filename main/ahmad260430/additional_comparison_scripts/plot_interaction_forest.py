#!/usr/bin/env python3
"""Forest plot of the factorial ANCOM-BC2 interaction coefficients.

The interaction terms are the only place the factorial fit answers "does
radiation change what the treatment does", and they are the one family of
coefficients that gets no `.res` file: a difference of differences has no
enriched class, so the cladogram/barplot encoding cannot carry it (see the
comment above the .res loop in run_ancombc2_factorial.R).

A forest plot can. Each taxon is one estimate with its 95% interval against a
zero line, so the width of the intervals -- the actual finding here -- is
visible rather than collapsed into a count of significant taxa.

The Day17-vs-Day25 joint fit (run_ancombc2_day17_day25_joint.R) has the same
kind of term, Timepoint:Arm, and is drawn with --dir/--prefix/--note.

Usage, inside the qiime2 container:
    python3 plot_interaction_forest.py [--db Greengenes2] [--suffix _Day17_minpres6]
    python3 plot_interaction_forest.py --dir comparison-day17-day25 \
        --prefix ANCOMBC2_day17day25_joint_ --suffix _armrefIFA \
        --note "Day17 vs Day25 · zero = both arms shift alike"
"""
import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.optimize import brentq
from scipy.stats import t as t_dist

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)

LN2 = math.log(2)
ALPHA = 0.05

# The study's PrimaryColor blues and oranges are tuned for filled areas -- as
# 6px dots and 1.5px intervals on white they fail contrast, lightness and
# chroma. These are darker steps of the same two hues that pass all six checks
# of the palette validator, so the treatment association survives at mark size.
TREATMENT_HUE = {"TreatmentIFA:RadiationIR": "#1F6FD0", "TreatmentInosine:RadiationIR": "#C25200",
                 "TimepointDay25:ArmIFA_IR": "#1F6FD0",
                 # The treatment-vs-treatment panel belongs to neither arm alone.
                 "IFA vs Inosine": "#4a3aa7"}
NULL_INK = "#9a9a97"
RULE_INK = "#c9c9c4"
TEXT_INK = "#2b2b28"
MUTED_INK = "#6b6b66"

# Named as the report table's interaction rows are, "A vs B": whether radiation
# moves the taxon differently in A than in B, with the first-named arm on the
# positive side.
PANEL_TITLE = {
    "TreatmentIFA:RadiationIR": "Isoferulic acid vs Control",
    "TreatmentInosine:RadiationIR": "Inosine vs Control",
    "TimepointDay25:ArmIFA_IR": "IFA+IR shift minus IFA shift",
    "IFA vs Inosine": "Isoferulic acid vs Inosine",
}

# The factorial forest shows the interaction for every pair of treatments: the
# two terms against Control from the main fit, and Inosine against IFA from the
# IFA-reference fit (read from `<stem>_refIFA.tsv`). The third panel adds no
# model information -- it is the IFA panel minus the Inosine panel -- but it can
# catch two arms moving in opposite directions, each too weakly to show alone.
# The fit's coefficient is Inosine relative to IFA; the sign of -1 turns it
# round so IFA is the first-named arm, as in the treatment rows.
THIRD_PANEL = ("IFA vs Inosine", "_refIFA", "TreatmentInosine:RadiationIR", -1)


def bh(p):
    """Benjamini-Hochberg q-values for a 1-D array of p-values (no NaNs)."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def joint_q(tables, panels):
    """q-values corrected across all panels at once, {panel: {taxon: q}}.

    The panels are one family, every pair of treatments, so the p-values of
    every taxon in every panel go through one Benjamini-Hochberg pass. Run on
    all tested taxa, before any genus-level filtering of what is drawn.
    """
    rows = [(panel, t, p) for panel, (fit, coef, _sign) in panels.items()
            for t, p in zip(tables[fit]["taxon"], tables[fit][f"p_{coef}"]) if pd.notna(p)]
    q = bh([r[2] for r in rows])
    out = {panel: {} for panel in panels}
    for (panel, t, _p), qv in zip(rows, q):
        out[panel][t] = qv
    return out


def load_panels(tab_dir, stem):
    """Tables keyed by fit suffix, and {panel: (fit suffix, coefficient, sign)}."""
    tsv = os.path.join(tab_dir, stem + ".tsv")
    if not os.path.exists(tsv):
        raise SystemExit(f"no such table: {tsv}")
    tables = {"": pd.read_csv(tsv, sep="\t")}
    panels = {c: ("", c, 1) for c in ("TreatmentIFA:RadiationIR", "TreatmentInosine:RadiationIR",
                                   "TimepointDay25:ArmIFA_IR") if f"lfc_{c}" in tables[""].columns}
    key, fit, coef, sign = THIRD_PANEL
    ref = os.path.join(tab_dir, f"{stem}{fit}.tsv")
    if coef in panels and os.path.exists(ref):
        tables[fit] = pd.read_csv(ref, sep="\t")
        panels[key] = (fit, coef, sign)
    return tables, panels


def short_name(lineage):
    """Deepest named rank, the way the barplots label a taxon."""
    parts = [p.strip() for p in str(lineage).split(";") if p.strip()]
    for p in reversed(parts):
        if not p.endswith("__"):
            return p.replace("_", " ")
    return parts[-1] if parts else str(lineage)


RANK_WORD = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
             "f": "family", "g": "genus", "s": "species"}


def axis_labels(lineages):
    """Deepest named rank per taxon, disambiguated where two lineages share one.

    A feature whose lineage stops at a genus and one that carries that genus
    plus an empty `s__` both reduce to the genus name, so without this the axis
    repeats a label and the two rows read as a duplicated taxon. The second is
    a real, distinct feature -- the genus with its species unresolved -- so say
    that rather than inventing a suffix.
    """
    names = [short_name(t) for t in lineages]
    clashes = {n for n in names if names.count(n) > 1}
    out = []
    for lineage, name in zip(lineages, names):
        parts = [p.strip() for p in str(lineage).split(";") if p.strip()]
        named = sum(1 for p in parts if not p.endswith("__"))
        if name not in clashes or named >= len(parts):
            out.append(name)
        else:
            rank = RANK_WORD.get(parts[named][:1], "rank")
            out.append(f"{name} · {rank} unresolved")
    return out


# Residual df of the factorial fit: 24 samples minus 6 coefficients. ANCOM-BC2's
# per-taxon df never exceeds it, so it is the fallback when W and p are too close
# to null to pin df down; such taxa cover zero at any df.
MAX_DF = 18


def implied_df(w, p):
    """The t degrees of freedom at which 2 * P(T > |w|) equals ANCOM-BC2's p."""
    if not (np.isfinite(w) and np.isfinite(p)) or p >= 0.999 or abs(w) < 0.05:
        return MAX_DF
    f = lambda v: 2 * t_dist.sf(abs(w), v) - p
    lo, hi = 0.05, 1e6
    if f(lo) * f(hi) > 0:
        return MAX_DF
    return brentq(f, lo, hi)


def coefficient_frame(df, coef, q_override=None, sign=1):
    """One row per taxon: estimate, interval, and which of the three states it is in."""
    lfc, se, q = f"lfc_{coef}", f"se_{coef}", f"q_{coef}"
    ss = f"passed_ss_{coef}"
    out = pd.DataFrame({
        "taxon": df["taxon"],
        "est": sign * pd.to_numeric(df[lfc], errors="coerce") / LN2,
        "se": pd.to_numeric(df[se], errors="coerce") / LN2,
        "W": sign * pd.to_numeric(df[f"W_{coef}"], errors="coerce"),
        "p": pd.to_numeric(df[f"p_{coef}"], errors="coerce"),
        "q": df["taxon"].map(q_override) if q_override is not None else pd.to_numeric(df[q], errors="coerce"),
        "ss": df[ss].astype(str).str.upper().eq("TRUE") if ss in df else True,
    })
    # ANCOM-BC2 tests W against a t distribution with per-taxon degrees of
    # freedom (1 to 18 here), not a normal, so +/-1.96 SE draws intervals that
    # clear zero for taxa whose own p is as high as 0.23. Recover each taxon's
    # df from its W and p and use the matching t quantile, so an interval
    # excludes zero exactly when ANCOM-BC2's p < 0.05.
    out["df"] = [implied_df(w, p) for w, p in zip(out.W, out.p)]
    half = t_dist.ppf(0.975, out.df) * out.se
    out["lo"] = out.est - half
    out["hi"] = out.est + half
    # Three states, in the order the reader should trust them.
    out["state"] = "null"
    out.loc[out.p < ALPHA, "state"] = "nominal"
    out.loc[(out.q < ALPHA) & out.ss, "state"] = "survives"
    return out.dropna(subset=["est", "se"])


def draw_clip_arrows(ax, frame, coef):
    """Mark intervals that run past the axis, so a clipped end is not read as the real bound."""
    x0, x1 = ax.get_xlim()
    for i, row in enumerate(frame.itertuples()):
        if not (np.isfinite(row.lo) and np.isfinite(row.hi)):
            continue
        ink = TREATMENT_HUE[coef] if row.state != "null" else NULL_INK
        if row.lo < x0:
            ax.plot([x0], [i], marker="<", ms=5, color=ink, clip_on=False, zorder=4)
        if row.hi > x1:
            ax.plot([x1], [i], marker=">", ms=5, color=ink, clip_on=False, zorder=4)


def draw_panel(ax, frame, coef, show_labels):
    hue = TREATMENT_HUE[coef]
    y = np.arange(len(frame))

    ax.axvline(0.0, color=TEXT_INK, lw=1.0, ls=(0, (4, 3)), zorder=1)

    for i, row in enumerate(frame.itertuples()):
        if not np.isfinite(row.est) or not np.isfinite(row.se):
            continue  # estimable in the other panel only; the row slot stays
        survives = row.state == "survives"
        nominal = row.state == "nominal"
        ink = hue if survives or nominal else NULL_INK
        ax.plot([row.lo, row.hi], [i, i], color=ink,
                lw=2.0 if survives else 1.4, solid_capstyle="round",
                alpha=1.0 if survives or nominal else 0.55, zorder=2)
        ax.plot([row.est], [i],
                marker="o", ms=7 if survives else 5.5,
                # Fill carries the state as well as colour does, so the three
                # groups stay apart in greyscale and for colour-blind readers.
                mfc=hue if survives else ("white" if nominal else NULL_INK),
                mec=ink, mew=1.6 if nominal else 1.0,
                alpha=1.0 if survives or nominal else 0.55, zorder=3)
        if survives:
            # Away from zero, so the label clears its own interval whichever
            # side of the zero line the estimate falls on.
            right = row.est > 0
            ax.annotate(f"q = {row.q:.1e}",
                        xy=(row.hi if right else row.lo, i),
                        xytext=(6 if right else -6, 0), textcoords="offset points",
                        va="center", ha="left" if right else "right",
                        fontsize=6.5, color=TEXT_INK)

    ax.set_yticks(y)
    if show_labels:
        ax.set_yticklabels(axis_labels(list(frame.taxon)), fontsize=6.5)
        for tick, state in zip(ax.get_yticklabels(), frame.state):
            tick.set_color(TEXT_INK if state == "survives" else MUTED_INK)
            if state == "survives":
                tick.set_fontweight("bold")
    else:
        ax.set_yticklabels([])
    ax.set_ylim(-1, len(frame))
    ax.invert_yaxis()

    ax.set_title(PANEL_TITLE[coef], fontsize=10, color=TEXT_INK, pad=10)
    ax.tick_params(axis="x", labelsize=7.5, colors=MUTED_INK)
    ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(RULE_INK)
    ax.grid(axis="x", color=RULE_INK, lw=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default="comparison-Treatment-x-Radiation-factorial")
    ap.add_argument("--prefix", default="ANCOMBC2_factorial_")
    ap.add_argument("--db", default="Greengenes2")
    ap.add_argument("--suffix", default="_Day17_minpres6")
    ap.add_argument("--note", default="Day17 · zero = radiation moves both treatments alike · "
                                      "right of zero = radiation raises the taxon more in the first-named")
    ap.add_argument("--genus-names-only", action="store_true",
                    help="draw only taxa whose deepest name is a genus or higher; same fit, "
                         "same q-values, species-resolved rows left out")
    args = ap.parse_args()

    tab_dir = os.path.join(STUDY_DIR, args.dir, "tables")
    fig_dir = os.path.join(STUDY_DIR, args.dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    stem = f"{args.prefix}{args.db}{args.suffix}"
    tables, panels = load_panels(tab_dir, stem)
    if not panels:
        raise SystemExit(f"{stem}.tsv carries no interaction coefficients")
    # Several panels are one family and are corrected together; a single panel
    # keeps the q-values as fitted.
    qj = joint_q(tables, panels) if len(panels) > 1 else None
    out_stem = stem
    if args.genus_names_only:
        # Not a re-test: q-values stay as corrected across every taxon; only rows
        # with a species name are left off the plot.
        for fit, d in tables.items():
            named = d["taxon"].map(lambda t: [p for p in str(t).split(";") if p.strip() and not p.strip().endswith("__")])
            tables[fit] = d[named.map(lambda ps: bool(ps) and not ps[-1].strip().startswith("s__"))].reset_index(drop=True)
        out_stem = f"{stem}_gplus"
    df = tables[""]
    coefs = list(panels)

    frames = {c: coefficient_frame(tables[panels[c][0]], panels[c][1], qj[c] if qj else None,
                                   sign=panels[c][2]).set_index("taxon")
              for c in coefs}
    # One shared row grid, so a row is the same taxon in every panel and the
    # reader scans across instead of re-reading names. Ordered by the first
    # term's estimate; a taxon estimable in only one panel keeps its slot in
    # the other and simply draws nothing there.
    union = sorted({t for f in frames.values() for t in f.index})
    order = list(frames[coefs[0]].est.reindex(union)
                 .sort_values(ascending=False, na_position="last").index)
    frames = {c: f.reindex(order).reset_index() for c, f in frames.items()}

    n = len(order)
    fig, axes = plt.subplots(
        1, len(coefs), sharex=True,
        figsize=(4.1 * len(coefs) + 2.6, max(4.0, 0.16 * n + 2.0)),
    )
    axes = np.atleast_1d(axes)
    # Header and footer sit a fixed distance in inches from the figure edge, so a
    # short figure (few taxa) does not stack them on top of each other; these
    # are the distances the original fractions gave at the 70-taxon height.
    H = fig.get_figheight()
    top_ = lambda inch: 1 - inch / H
    bot_ = lambda inch: inch / H

    for ax, coef in zip(axes, coefs):
        draw_panel(ax, frames[coef], coef, show_labels=(ax is axes[0]))

    # Panels share one x scale: the interaction terms are the same quantity
    # and must be comparable by eye.
    # A taxon with df near 1 has a t quantile near 12.7 and an interval tens of
    # log2 units wide; scaling to it would flatten every other row. Scale to the
    # taxa with df >= 3 and let draw_panel clip the rest with an arrow.
    lo = min(f.lo[f.df >= 3].min(skipna=True) for f in frames.values())
    hi = max(f.hi[f.df >= 3].max(skipna=True) for f in frames.values())
    pad = 0.08 * (hi - lo)
    for ax in axes:
        ax.set_xlim(lo - pad, hi + pad)
    for ax, coef in zip(axes, coefs):
        draw_clip_arrows(ax, frames[coef], coef)

    fig.suptitle("ANCOM-BC2 interaction coefficients", fontsize=12.5,
                 color=TEXT_INK, y=top_(0.20))
    note_head, _, note_tail = args.note.partition(" · ")
    fig.text(0.5, top_(0.59), f"{args.db} · {note_head} · {'genus-level names only · ' if args.genus_names_only else ''}"
             f"95% CI from ANCOM-BC2's own t test · {'q corrected across all panels · ' if qj else ''}{note_tail}", ha="center", fontsize=8, color=MUTED_INK)
    fig.text(0.5, bot_(0.82), "Interaction log2 fold change (95% CI)",
             ha="center", fontsize=9, color=TEXT_INK)

    survives = sum(int((f.state == "survives").sum()) for f in frames.values())
    nominal = sum(int((f.state == "nominal").sum()) for f in frames.values())
    null = sum(int((f.state == "null").sum()) for f in frames.values())
    handles = [
        Line2D([], [], marker="o", ls="none", ms=7, mfc=TEXT_INK, mec=TEXT_INK,
               label=f"survives BH correction, q < {ALPHA} ({survives})"),
        Line2D([], [], marker="o", ls="none", ms=5.5, mfc="white", mec=TEXT_INK,
               mew=1.6, label=f"interval excludes zero, not a corrected hit ({nominal})"),
        Line2D([], [], marker="o", ls="none", ms=5.5, mfc=NULL_INK, mec=NULL_INK,
               label=f"interval covers zero ({null})"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, bot_(0.40)),
               labelcolor=MUTED_INK, handletextpad=0.5, columnspacing=1.8)

    dropped = len(df) - len(order)
    if dropped > 0:
        shown = "genus-level-named taxa" if args.genus_names_only else "tested taxa"
        fig.text(0.5, bot_(0.16), f"{dropped} of {len(df)} {shown} have no estimable "
                 "standard error in any panel and are not shown.",
                 ha="center", fontsize=7, color=MUTED_INK)

    fig.subplots_adjust(top=top_(0.92), bottom=bot_(1.25), left=0.30, right=0.97, wspace=0.10)

    for ext in ("png", "svg"):
        out = os.path.join(fig_dir, f"{out_stem}_interaction_forest.{ext}")
        fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.15)
        print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)
    print(f"  {n} taxa plotted · {survives} survive · {nominal} nominal · {null} null")


if __name__ == "__main__":
    main()
