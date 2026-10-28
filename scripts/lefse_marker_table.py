"""Write a LEfSe results file's markers as a readable table.

One row per marker (a taxon LEfSe assigned to a group), in the given group
order and then by LDA score, largest first. Taxa that were tested but not
assigned stay in the raw results file only.

Run by the lefse_Huttenhower_table rule; project scripts import
write_lefse_marker_table() directly.
"""

import csv
import re

RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
# LEfSe's input step names an unassigned rank after the deepest real call plus
# the tree level it stands in for: "g__Duncaniella_s__L7", "f__Muribaculaceae_x__L6".
PLACEHOLDER = re.compile(r"^(.+)_[a-z]__L(\d)$")

COLUMNS = ["group", "taxon", "tree_level", "lda_score",
           "kruskal_wallis_p_uncorrected", "log10_highest_group_mean", "lineage"]


def write_lefse_marker_table(res_path, out_path, group_order):
    rows = []
    with open(res_path) as fh:
        for line in fh:
            lineage, logmax, group, lda, p = line.rstrip("\n").split("\t")
            if not group:
                continue
            if group not in group_order:
                raise ValueError(f"Group {group!r} in {res_path} is not in {group_order}")
            ranks = lineage.split(".")
            name = ranks[-1]
            m = PLACEHOLDER.match(name)
            if m:
                name = f"{m.group(1)} (no {RANKS[int(m.group(2)) - 1].lower()} call)"
            rows.append([group, name, RANKS[len(ranks) - 1], f"{float(lda):.3f}", p,
                         f"{float(logmax):.3f}", "; ".join(ranks)])
    rows.sort(key=lambda r: (group_order.index(r[0]), -float(r[3])))
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)
    return len(rows)


if "snakemake" in globals():
    n = write_lefse_marker_table(snakemake.input.res, snakemake.output.table,
                                 list(snakemake.params.group_order))
    with open(snakemake.log[0], "w") as log:
        log.write(f"Wrote {n} markers to {snakemake.output.table}\n")
