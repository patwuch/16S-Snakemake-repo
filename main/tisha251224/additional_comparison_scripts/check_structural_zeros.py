"""Report check for the structural-zero presence table, shared by the report checkers.

check(html, text) compares the presence table (the one headed "Absent from first
/ second") row by row with structural_zeros.summary() and checks the headline
numbers the prose quotes. Structural zeros are described, not tested (2026-10-01):
the table carries no p or q, and the report quotes none. check_dots(html, out_dir)
checks the dot map's caption against its sidecar. Each returns a list of misses.
"""
import os
import re

import pandas as pd

import structural_zeros as sz


def check(html, text):
    misses = []
    df = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros.tsv"), sep="\t")
    fresh = sz.tests()
    if len(fresh) != len(df) or not (fresh.p.round(12).values == df.p.round(12).values).all():
        misses.append("structural_zeros.tsv is stale: rerun structural_zeros.py")
    for c in ("carrier_median_pct", "carrier_min_pct", "carrier_max_pct"):
        if not (fresh[c].round(12).values == df[c].round(12).values).all():
            misses.append(f"structural_zeros.tsv {c} is stale: rerun structural_zeros.py")
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
    if re.search(r"Fisher|Firth|calibrat|binary.Jaccard|presence (test|interaction|call)|permutation estimate", text):
        misses.append("the report still cites a presence statistic; structural zeros are described, not tested")
    # Carriers' share of reads: the numbers the prose quotes.
    for w in (f"The median over all {len(df)} calls is {sz.pct_str(df.carrier_median_pct.median())}% of a sample's reads, "
              f"{int((df.carrier_median_pct < 0.1).sum())} are below 0.1%, and in {int((df.carrier_max_pct < 0.1).sum())} "
              f"even the richest carrier stays below 0.1%",
              f"One read is at most {sz.pct_str(sz.one_read_pct())}% of any sample's reads"):
        if w not in text:
            misses.append(f"carrier share: '{w}' not in report")
    for t, g, lead, tail in (("Amulumruptor caecigallinarius", "RS+IR", "holds ", " of the reads in the three RS+IR samples"),
                             ("Mucispirillum schaedleri", "Cellulose", "holds ", " in the three Cellulose samples")):
        r = df[df.taxon.str.endswith(t) & (df.present_in == g)].iloc[0]
        w = f"{lead}{r.carrier_share}{tail}"
        if w not in text or r.n_present != 3:
            misses.append(f"carrier share {t} in {g}: '{w}' not in report")
    # Complete separations: the taxon in every sample of one group and none of the other.
    cs = df[df.n_present == df.n_group]
    n = lambda c: len(cs[cs.contrast == c])
    m = lambda c, g: int((cs[cs.contrast == c].absent_from == g).sum())
    for w in (f"{len(cs)} calls, involving {cs.taxon.nunique()} taxa, are complete separations",
              f"{n('RS+IR vs Cellulose+IR')} separate RS+IR from Cellulose+IR, {m('RS+IR vs Cellulose+IR', 'RS+IR')} of them taxa missing from RS+IR",
              f"{n('Inulin+IR vs Cellulose+IR')} separate Inulin+IR from Cellulose+IR, {m('Inulin+IR vs Cellulose+IR', 'Inulin+IR')} of them missing from Inulin+IR",
              f"{n('RS+IR vs RS')} separate RS+IR from unirradiated RS, {m('RS+IR vs RS', 'RS+IR')} of them missing after irradiation",
              f"inulin and cellulose separate on {n('Inulin vs Cellulose')} taxa and RS and cellulose on {n('RS vs Cellulose')}",
              f"{n('Inulin vs Cellulose')} taxa separate inulin from cellulose completely and {n('RS vs Cellulose')} RS from",
              f"{n('RS+IR vs Cellulose+IR')} taxa separate RS+IR from Cellulose+IR completely, {m('RS+IR vs Cellulose+IR', 'RS+IR')} of them",
              f"and {n('Inulin+IR vs Cellulose+IR')} separate Inulin+IR, {m('Inulin+IR vs Cellulose+IR', 'Inulin+IR')} of them missing from Inulin+IR",
              f"{m('RS+IR vs RS', 'RS+IR')} of the {n('RS+IR vs RS')} taxa that separate RS+IR from unirradiated RS completely"):
        if w not in text:
            misses.append(f"complete separations: '{w}' not in report")
    return misses, df


def check_dots(html, out_dir):
    """The dot map's caption quotes its sidecar (plot_presence_dots_factorial.py)."""
    misses = []
    for gplus in ("", "_gplus"):
        stem = f"{sz.STEM}{gplus}_presence_dots"
        m = re.search(rf'src="figures/{stem}\.png"', html)
        if not m:
            misses.append(f"report does not show {stem}")
            continue
        s = pd.read_csv(os.path.join(out_dir, "figures", f"{stem}.tsv"), sep="\t").iloc[0]
        s = {k: (int(v) if k not in ("min_pct", "max_pct") else v) for k, v in s.items()}
        cap = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html[m.end():html.index("</figcaption>", m.end())]))
        if gplus:
            want = [f"for the {s['taxa']} taxa named no further than genus; {s['structural_zeros']} of them are structural zeros"]
        else:
            ab = {g: int(s[f"absent_{g}"]) for g in ("RS+IR", "RS", "Inulin+IR", "Inulin", "Cellulose", "Cellulose+IR")}
            if ab["Cellulose"] != ab["Cellulose+IR"]:
                misses.append(f"{stem}: the caption gives one count for both cellulose groups, the sidecar {ab}")
            want = [f"the {s['taxa'] - s['structural_zeros']} seen in every group and the {s['structural_zeros']} structural zeros: {s['taxa']} in all",
                    f"{ab['RS+IR']} in RS+IR, {ab['RS']} in RS, {ab['Inulin+IR']} in Inulin+IR, {ab['Inulin']} in Inulin "
                    f"and {ab['Cellulose']} each in Cellulose and Cellulose+IR"]
        for w in want:
            if w not in cap:
                misses.append(f"{stem}: caption lacks '{w}' from its sidecar")
    return misses


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
