"""Report check for the structural-zero presence table, shared by the report checkers.

check(html, text[, tests]) compares the presence table (the one headed "Absent from first
/ second") row by row with structural_zeros.summary() and checks the headline
numbers the prose quotes. Returns a list of misses.
"""
import os
import re

import pandas as pd

import structural_zeros as sz


def _fmt(v):
    return "—" if v is None or v != v else (f"{v:.2g}" if v >= 0.01 else f"{v:.1e}")


def check(html, text, tests=True):
    """tests=False: the table carries only the counts (no min p / min q) and the prose quotes no test."""
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
        if len(c) == (7 if tests else 5):
            question, c = c[0], c[1:]
        got[(question, c[0])] = tuple(c[1:])
    for r in sz.summary(df):
        want = (str(r["n"]), f"{r['absent_first']} / {r['absent_second']}", str(r["n_half"]),
                _fmt(r["min_p"]), _fmt(r["min_q"]))[:None if tests else 3]
        if got.get((r["question"], r["label"])) != want:
            misses.append(f"presence row {r['label']}: report {got.get((r['question'], r['label']))}, tables {want}")
    if not tests:
        return misses, df
    if f"across all {len(df)} such tests" not in text:
        misses.append(f"presence test count {len(df)} not in report")
    if (df.q < 0.05).any():
        misses.append("a structural zero has q < 0.05 but the report says none does")
    if f"the smallest q is {df.q.min():.2f}" not in text:
        misses.append(f"smallest presence q {df.q.min():.2f} not in report")
    return misses, df


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
