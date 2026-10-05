#!/usr/bin/env python3
"""Refresh every file built from a report.html, in order, after any report change.

  1. the report's checker (qiime2 container): stop here if any quoted number is wrong
  2. report_figures.zip       make_report_figures_zip.py, for the report that offers it
  3. report-standalone.html   report.html with every figure inlined as a base64 PNG, one file
                              that opens anywhere; any "Download all figures" row is dropped,
                              since its link cannot resolve from a single file
  4. report-share.zip         report.html, report_figures.zip where offered, and the PNGs it
                              shows, laid out so the report opens with working links once unzipped
  5. report.pdf               report_to_pdf.py (Playwright): every panel view shown and labelled

Each output is written to a temporary name and moved into place only when complete,
so a failed step never leaves a half-written file behind.

Usage, with the host python3 (Playwright for the PDF; apptainer for the checker):
    python3 build_report_outputs.py [factorial] [day17-day25] [--skip-check] [--skip-pdf]
With no report named, both are rebuilt.
"""
import argparse
import base64
import os
import re
import subprocess
import sys
import zipfile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
REPO = os.path.dirname(os.path.dirname(STUDY_DIR))
QIIME = os.path.join(REPO, "containers", "qiime2-2025.10-amplicon-core.sif")
DL_ROW = re.compile(r'\s*<p class="dl-row">.*?</p>', re.S)
REPORTS = {
    "factorial": dict(dir="comparison-Treatment-x-Radiation-factorial", checker="check_report_factorial.py",
                      figure_zip=True),
    "day17-day25": dict(dir="comparison-day17-day25", checker="check_report_day17_day25.py", figure_zip=False),
}


def step(name):
    print(f"== {name}", flush=True)


def figures(html):
    return sorted(set(re.findall(r'src="(figures/[^"]+\.png)"', html)))


def replace_atomically(tmp, final):
    os.replace(tmp, final)
    print(f"   wrote {os.path.relpath(final, STUDY_DIR)} ({os.path.getsize(final) / 1e6:.1f} MB)")


def run_check(checker):
    step(checker)
    r = subprocess.run(["apptainer", "exec", "--bind", f"{REPO}:{REPO}", "--pwd", THIS_DIR, QIIME,
                        "python3", checker], capture_output=True, text=True)
    out = "\n".join(l for l in (r.stdout + r.stderr).splitlines() if not l.startswith("INFO:"))
    print("   " + out.replace("\n", "\n   "))
    if r.returncode != 0:
        sys.exit(f"{checker} failed: nothing rebuilt")


def build_standalone(out_dir, html):
    step("report-standalone.html")
    figs = figures(html)
    missing = [f for f in figs if not os.path.exists(os.path.join(out_dir, f))]
    if missing:
        sys.exit(f"missing figures: {missing[:3]}")
    out = DL_ROW.sub("", html, count=1)
    for f in figs:
        data = base64.b64encode(open(os.path.join(out_dir, f), "rb").read()).decode("ascii")
        out = out.replace(f'src="{f}"', f'src="data:image/png;base64,{data}"')
    left = re.findall(r'src="figures/', out)
    if left:
        sys.exit(f"{len(left)} figure reference(s) not inlined")
    if 'href="report_figures.zip"' in out:
        sys.exit("download link survived in the standalone report")
    tmp = os.path.join(out_dir, ".report-standalone.html.tmp")
    open(tmp, "w").write(out)
    replace_atomically(tmp, os.path.join(out_dir, "report-standalone.html"))
    print(f"   {len(figs)} figures inlined")


def build_share(out_dir, html, figure_zip):
    step("report-share.zip")
    figs = figures(html)
    tmp = os.path.join(out_dir, ".report-share.zip.tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(os.path.join(out_dir, "report.html"), "report.html")
        if figure_zip:
            z.write(os.path.join(out_dir, "report_figures.zip"), "report_figures.zip")
        for f in figs:
            z.write(os.path.join(out_dir, f), f)
    replace_atomically(tmp, os.path.join(out_dir, "report-share.zip"))
    print(f"   report.html{', report_figures.zip' if figure_zip else ''} and {len(figs)} figures")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="*", help=f"reports to rebuild: {', '.join(REPORTS)} (default: both)")
    ap.add_argument("--skip-check", action="store_true", help="do not run the report checker first")
    ap.add_argument("--skip-pdf", action="store_true", help="do not reprint report.pdf")
    args = ap.parse_args()
    names = args.reports or list(REPORTS)
    unknown = [n for n in names if n not in REPORTS]
    if unknown:
        ap.error(f"unknown report(s): {', '.join(unknown)}")

    if not args.skip_check:
        for name in names:
            run_check(REPORTS[name]["checker"])
    for name in names:
        spec = REPORTS[name]
        out_dir = os.path.join(STUDY_DIR, spec["dir"])
        print(f"### {spec['dir']}")
        if spec["figure_zip"]:
            step("report_figures.zip")
            subprocess.run([sys.executable, os.path.join(THIS_DIR, "make_report_figures_zip.py")], check=True)
        html = open(os.path.join(out_dir, "report.html")).read()
        build_standalone(out_dir, html)
        build_share(out_dir, html, spec["figure_zip"])
        if not args.skip_pdf:
            step("report.pdf")
            subprocess.run([sys.executable, os.path.join(THIS_DIR, "report_to_pdf.py"),
                            os.path.join(out_dir, "report.html"), os.path.join(out_dir, "report.pdf")], check=True)
    print("all outputs refreshed")


if __name__ == "__main__":
    main()
