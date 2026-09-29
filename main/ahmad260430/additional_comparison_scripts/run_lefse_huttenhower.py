#!/usr/bin/env python3
"""Run the Huttenhower LEfSe implementation on a pre-filtered feature table.

Standalone counterpart to the shared pipeline's scripts/run_lefse_Huttenhower.py
(which takes its inputs from a Snakemake `script:` binding). Same three-step
process (QIIME2 artifacts -> LEfSe TSV -> .in -> .res), just argparse'd instead.
Must run inside the LEfSe container (needs its own python + R env, and
scripts/lefse/*.py from the main repo).
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

COLLAPSE_LEVEL = {
    "Silva138":    6,   # genus
    "Greengenes2": 7,   # species
}


def parse_colors(pairs):
    colors = {}
    for pair in pairs:
        name, hexval = pair.split("=", 1)
        colors[name] = hexval
    return colors


def run(cmd):
    print("Running:", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-table", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--group-col", required=True)
    parser.add_argument("--db", required=True, choices=sorted(COLLAPSE_LEVEL))
    parser.add_argument("--lda-cutoff", type=float, required=True)
    parser.add_argument("--kw-cutoff", type=float, required=True)
    parser.add_argument("--wilcoxon-cutoff", type=float, required=True)
    parser.add_argument("--color", action="append", required=True, metavar="GROUP=HEX",
                         help="repeatable; order given fixes the LDA/cladogram group order")
    parser.add_argument("--lefse-dir", required=True, help="path to scripts/lefse/")
    parser.add_argument("--output", required=True, help="path for the .res/.txt result file")
    args = parser.parse_args()

    colors_map = parse_colors(args.color)
    collapse_level = COLLAPSE_LEVEL[args.db]
    lefse_dir = Path(args.lefse_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def lefse_script(name):
        p = lefse_dir / name
        if not p.exists():
            raise FileNotFoundError(f"LEfSe script not found: {p}")
        return str(p)

    order = list(colors_map.keys())
    sortkey = {}
    for i, g in enumerate(order):
        g_norm = g.replace(" ", "_")
        label = f"zzsort{i:02d}_{g_norm}"
        sortkey[g] = label
        sortkey[g_norm] = label

    def strip_prefix(path):
        p = Path(path)
        text = p.read_text()
        for g in order:
            text = text.replace(sortkey[g], g)
        p.write_text(text)

    with tempfile.TemporaryDirectory() as _tmp:
        tmpdir = Path(_tmp)
        lefse_input_tsv = tmpdir / "lefse_input.tsv"
        lefse_input_in = tmpdir / "lefse_input.in"

        run([
            sys.executable, lefse_script("plugin_lefse_input.py"),
            "--table_file", args.feature_table,
            "--taxonomy_file", args.taxonomy,
            "--metadata_file", args.metadata,
            "--output_file", str(lefse_input_tsv),
            "--class_col", args.group_col,
            "--level", str(collapse_level),
        ])

        tsv = pd.read_csv(lefse_input_tsv, sep="\t", header=None,
                           dtype=str, keep_default_na=False)
        feature_col = tsv.iloc[:, 0]
        data = tsv.iloc[:, 1:].copy()

        present = set(data.iloc[0].unique())
        missing = present - set(sortkey)
        if missing:
            raise ValueError(
                f"Group(s) in the data are absent from --color: {sorted(missing)}.\n"
                f"--color keys given: {sorted(sortkey)}"
            )

        data.iloc[0] = data.iloc[0].map(sortkey)
        pd.concat([feature_col, data], axis=1).to_csv(
            lefse_input_tsv, sep="\t", header=False, index=False
        )

        run([
            sys.executable, lefse_script("plugin_lefse_format.py"),
            str(lefse_input_tsv), str(lefse_input_in),
            "-c", "1", "-o", "1000000",
        ])

        run([
            sys.executable, lefse_script("plugin_lefse_run.py"),
            str(lefse_input_in), str(output_path),
            "-l", str(args.lda_cutoff),
            "-a", str(args.kw_cutoff),
            "-w", str(args.wilcoxon_cutoff),
        ])

        strip_prefix(output_path)

    print(f"Wrote LEfSe results -> {output_path}")


if __name__ == "__main__":
    main()
