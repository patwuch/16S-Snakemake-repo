"""Report check for the structural-zero presence table, shared by the report checkers.

check(html, text) compares the presence table (the one headed "Absent from first
/ second") row by row with structural_zeros.summary() and checks the headline
numbers the prose quotes. Structural zeros are described, not tested (2026-10-01):
the table carries no p or q, and the report quotes none. Returns a list of misses.
"""
import os
import re

import pandas as pd

import structural_zeros as sz


def _fmt(v):
    return "—" if v is None or v != v else (f"{v:.2g}" if v >= 0.01 else f"{v:.1e}")


def check(html, text):
    misses = []
    df = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros.tsv"), sep="\t")
    fresh = sz.tests()
    if len(fresh) != len(df) or not (fresh.p.round(12).values == df.p.round(12).values).all():
        misses.append("structural_zeros.tsv is stale: rerun structural_zeros.py")
    i = html.index("Absent from first / second")
    body = html[html.index("<tbody>", i):html.index("</tbody>", i)]
    got, question = {}, None
    for r in re.findall(r"<tr[^>]*>(.*?)</tr>", body, re.S):
        c = [re.sub(r"<[^>]+>", "", x).strip() for x in re.findall(r"<td[^>]*>(.*?)</td>", r)]
        if len(c) == 6:
            question, c = c[0], c[1:]
        got[(question, c[0])] = tuple(c[1:])
    for r in sz.summary(df):
        want = (str(r["n"]), f"{r['absent_first']} / {r['absent_second']}", str(r["n_half"]),
                f"{sz.pct_str(r['med_lo'])}–{sz.pct_str(r['med_hi'])}%")
        if got.get((r["question"], r["label"])) != want:
            misses.append(f"presence row {r['label']}: report {got.get((r['question'], r['label']))}, tables {want}")
    if re.search(r"Fisher|Firth|calibrat|binary.Jaccard|presence (test|interaction|call)", text):
        misses.append("the report still cites a presence statistic; structural zeros are described, not tested")
    # Carriers' share of reads: recomputed from the counts, and the numbers the prose quotes.
    for c in ("carrier_median_pct", "carrier_min_pct", "carrier_max_pct"):
        if not (fresh[c].round(12).values == df[c].round(12).values).all():
            misses.append(f"structural_zeros.tsv {c} is stale: rerun structural_zeros.py")
    for w in (f"The median over all {len(df)} calls is {sz.pct_str(df.carrier_median_pct.median())}% of a mouse's reads, "
              f"{int((df.carrier_median_pct < 0.1).sum())} are below 0.1%, and in {int((df.carrier_max_pct < 0.1).sum())} "
              f"even the richest carrier stays below 0.1%",
              f"One read is at most {sz.pct_str(sz.one_read_pct())}% of any mouse's reads"):
        if w not in text:
            misses.append(f"carrier share: '{w}' not in report")
    for t, g, lead, tail in (("s__Bacteroides_H_857956 rodentium", "SCFAs", "Bacteroides_H rodentium at ", " in the {n} unirradiated SCFAs mice"),
                             ("s__Bacteroides_H_857956 rodentium", "IR", "(", " in {n} IR mice)"),
                             ("s__Ruminococcus_C_59129 flavefaciens_B", "Control", "flavefaciens_B at ", " in {n} control-diet mice"),
                             ("s__Ruminococcus_C_59129 flavefaciens_B", "Control", "4 of 8 control-diet mice (", " of their reads)"),
                             ("s__Ruminococcus_C_59129 flavefaciens_B", "SCFAs+IR", "4 of 8 irradiated SCFAs mice (", ")"),
                             ("g__COE1;s__", "SCFAs", "3 of 8 unirradiated SCFAs mice (", ")"),
                             ("g__COE1;s__", "IR", "5 of 8 irradiated control-diet mice (", ")")):
        r = df[df.taxon.str.endswith(t) & (df.present_in == g)].iloc[0]
        w = f"{lead}{r.carrier_share}{tail.format(n=r.n_present)}"
        if w not in text:
            misses.append(f"carrier share {t} in {g}: '{w}' not in report")
    return misses, df


def check_heattrees(html, out_dir):
    """Every heat-tree view quotes its sidecar's counts (plot_heat_trees_factorial.R).
    Returns (misses, n_views)."""
    misses, n = [], 0
    for m in re.finditer(r'src="figures/([^"]+_heattree_[^"]+)\.png"', html):
        stem = m.group(1)
        side = os.path.join(out_dir, "figures", f"{stem}.tsv")
        if not os.path.exists(side):
            misses.append(f"{stem}: heat-tree sidecar missing")
            continue
        rows = [l.rstrip("\n").split("\t") for l in open(side)]
        s = dict(zip(rows[0], rows[1]))
        # The key list beside the tree must hold exactly the ringed taxa.
        keyed = s["q_hits"] if "q_hits" in s else s["labelled_hits"]
        # Panels ring their structural zeros; the interaction tree rings every untested tip.
        ringed = int(keyed) + int(s.get("no_estimate", s.get("structural_zeros", 0)))
        kpath = os.path.join(out_dir, "figures", f"{stem}.keys.tsv")
        nkeys = sum(1 for _ in open(kpath)) - 1 if os.path.exists(kpath) else 0
        if nkeys != ringed:
            misses.append(f"{stem}: key list has {nkeys} entries, expected {ringed} ringed taxa")
        # A structural zero's key entry ends with its presence counts and carriers' share, as tabled.
        if not stem.endswith("_interaction") and os.path.exists(kpath):
            got = sorted(l.rstrip("\n").split("\t")[2].rsplit("   ", 1)[1] for l in list(open(kpath))[1:]
                         if l.split("\t")[1] == "sz")
            szt = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros.tsv"), sep="\t")
            szt = szt[szt.contrast == f"{s['pos']} vs {s['neg']}"]
            want = sorted((f"{r.n_present}/{r.n_group} vs 0/{r.n_group}" if r.present_in == s["pos"]
                           else f"0/{r.n_group} vs {r.n_present}/{r.n_group}") + f" · {r.carrier_share}"
                          for r in szt.itertuples())
            # The genus view keeps a subset of the taxa, so its entries need only be among them.
            if ("_gplus_" not in stem and got != want) or ("_gplus_" in stem and not set(got) <= set(want)):
                misses.append(f"{stem}: structural-zero key entries differ from structural_zeros.tsv")
        cap = html[m.end():html.index("</figcaption>", m.end())]
        cap = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", cap))
        if stem.endswith("_interaction"):
            # Each dotted ring is keyed with its carriers per group, as the untestable table has them.
            un = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_untestable_absent_from_a_cell.tsv"), sep="\t")
            cells = [c for c in un.columns if " (n=" in c]
            carriers = {" · ".join(f"{c.split(' (n=')[0]} {v}/{c.split('n=')[1].rstrip(')')}" for c, v in zip(cells, r))
                        for r in un[cells].itertuples(index=False)}
            got = [l.rstrip("\n").split("\t")[2].rsplit("   ", 1)[1] for l in list(open(kpath))[1:] if l.split("\t")[1] == "sz"]
            if len(got) != int(s["no_estimate"]) or not set(got) <= carriers:
                misses.append(f"{stem}: dotted-ring key entries differ from the untestable table")
            if "_gplus_" in stem:
                want = [f"limited to the {s['tips']} tips", f"{s['estimated']} with an abundance estimate",
                        f"{s['q_hits']} at q &lt; 0.05", f"{s['no_estimate']} structural zeros with their carriers"]
            else:
                # Tip size: mice carrying the taxon. The caption's ranges, recounted from the tables.
                ix = pd.read_csv(os.path.join(sz.TAB_DIR, "DA_interaction_intersection_Greengenes2.tsv"), sep="\t")
                hits = 32 - ix.n_zero_mice[ix.ancombc2_q < 0.05]
                part = hits[hits < 32]
                want = [f"Each of the {s['estimated']} taxa seen in every group", f"The {s['no_estimate']} taxa absent",
                        f"the {s['q_hits']} taxa at q &lt; 0.05", f"{s['tips']} tips in all",
                        f"from {s['min_mice']} to {s['max_mice']} (key at the bottom right)",
                        f"(32 mice) are carried by only {part.min()} to {part.max()}",
                        f"and the structural zeros by {un.n_samples.min()} to {un.n_samples.max()}"]
                if int(s["hits_in_all_mice"]) != 1 or (hits == 32).sum() != 1:
                    misses.append(f"{stem}: caption says one hit is in all 32 mice; the tables give {(hits == 32).sum()}")
            for w in want:
                if w not in cap:
                    misses.append(f"{stem}: caption lacks '{w}' from its sidecar")
            n += 1
            continue
        where = " and ".join(p for p in (f"{s['sz_in_pos']} in {s['pos']}" if int(s["sz_in_pos"]) else "",
                                         f"{s['sz_in_neg']} in {s['neg']}" if int(s["sz_in_neg"]) else "") if p)
        noun = "taxon" if int(s["structural_zeros"]) == 1 else "taxa"
        if f"Also drawn: {s['structural_zeros']} structural-zero {noun}, seen in one group only ({where})" not in cap:
            misses.append(f"{stem}: caption structural zeros do not match its sidecar ({s['structural_zeros']}; {where})")
        if "_gplus_" in stem:
            want = [f"limited to the {s['tips']} tips", f"{s['labelled_hits']} reach q &lt; 0.05"]
        else:
            want = [f"{s['labelled_hits']} taxa reach q &lt; 0.05", f"{s['tips']} tips in all"]
        for w in want:
            if w not in cap:
                misses.append(f"{stem}: caption lacks '{w}' from its sidecar")
        n += 1
    return misses, n


def check_figures(html, out_dir):
    """Every view drawn from an _sz figure quotes its sidecar's counts, and no panel
    that has an _sz figure still shows the tested-only one. Returns (misses, n_views)."""
    misses, n = [], 0
    for m in re.finditer(r'src="figures/([^"]+)_cladogram\.png"', html):
        stem = m.group(1)
        side = os.path.join(out_dir, "tables", f"{stem}.struczero.tsv")
        if "_sz" not in stem:
            cands = [c for c in os.listdir(os.path.join(out_dir, "tables"))
                     if c.endswith(".struczero.tsv") and c.replace("_sz", "", 1).replace(".struczero.tsv", "") == stem]
            if cands:
                misses.append(f"{stem}: an _sz figure exists but the report shows the tested-only one")
            continue
        rows = [l.rstrip("\n").split("\t") for l in open(side)]
        a, b = rows[0][1], rows[0][2]
        rows = rows[1:]
        na, nb = sum(r[1] == a for r in rows), sum(r[1] == b for r in rows)
        where = " and ".join(p for p in (f"{na} in {a}" if na else "", f"{nb} in {b}" if nb else "") if p)
        cap = html[m.end():html.index("</figcaption>", m.end())]
        cap = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", cap))
        noun = "taxon" if len(rows) == 1 else "taxa"
        if f"Also drawn: {len(rows)} structural-zero {noun}, seen in one group only ({where})" not in cap:
            misses.append(f"{stem}: caption does not match its sidecar ({len(rows)}; {where})")
        n += 1
    return misses, n
