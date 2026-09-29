#!/usr/bin/env python3
"""Pack every figure report.html shows into report_figures.zip beside it.

The report's download button links to this file. The figure list is read from
the report's own <img> tags, so the zip always matches what the page shows;
each PNG goes in with its SVG twin where one exists. Re-run after redrawing
figures or changing which ones the report uses.

Usage:
    python3 make_report_figures_zip.py
"""
import os
import re
import zipfile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
REPORT = os.path.join(OUT_DIR, "report.html")
ZIP = os.path.join(OUT_DIR, "report_figures.zip")


def main():
    html = open(REPORT).read()
    pngs = sorted(set(re.findall(r'src="(figures/[^"]+\.png)"', html)))
    files = []
    for png in pngs:
        files.append(png)
        svg = png[:-4] + ".svg"
        if os.path.exists(os.path.join(OUT_DIR, svg)):
            files.append(svg)
    missing = [f for f in files if not os.path.exists(os.path.join(OUT_DIR, f))]
    if missing:
        raise SystemExit(f"report references figures that do not exist: {missing}")
    tmp = ZIP + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for f in files:
            z.write(os.path.join(OUT_DIR, f), arcname=os.path.join("report_figures", os.path.basename(f)))
    os.replace(tmp, ZIP)
    print(f"wrote {ZIP}: {len(pngs)} figures, {len(files)} files, {os.path.getsize(ZIP) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
