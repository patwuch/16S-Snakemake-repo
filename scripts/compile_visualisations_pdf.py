import os
# Send this job's stdout/stderr, subprocesses included, to its Snakemake log.
_log = os.open(snakemake.log[0], os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
os.dup2(_log, 1)
os.dup2(_log, 2)

import csv
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import cairosvg
import pypdf

# ── Snakemake inputs / params ─────────────────────────────────────
plot_inputs    = [Path(p) for p in snakemake.input.plots]
metadata_path  = Path(snakemake.input.metadata)
permanova_paths = [Path(p) for p in snakemake.input.permanova]
rarefy_base    = Path(snakemake.input.rarefy_base)
rarefy_dbs     = [Path(p) for p in snakemake.input.rarefy_dbs]
reference_dbs  = list(snakemake.params.reference_dbs)
phylogenetic_reference_dbs = list(snakemake.params.phylogenetic_reference_dbs)
trunc_len_csv  = snakemake.params.get("trunc_len_csv")   # None for TGS
out_pdf        = Path(str(snakemake.output.pdf))

A4_W, A4_H = 11.69, 8.27   # landscape A4, inches

HEADER_BG   = "#2c3e50"
HEADER_FG   = "white"
ALT_ROW     = "#f5f6fa"
ACCENT      = "#2980b9"
VALUE_COLOR = "#16a085"

# ── Core helpers ──────────────────────────────────────────────────

def add_figure(writer: pypdf.PdfWriter, fig: plt.Figure) -> None:
    """Convert a matplotlib figure to a PDF page and append it."""
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    for page in pypdf.PdfReader(buf).pages:
        writer.add_page(page)


def add_pdf_file(writer: pypdf.PdfWriter, path: Path) -> None:
    try:
        for page in pypdf.PdfReader(str(path)).pages:
            writer.add_page(page)
    except Exception as exc:
        print(f"[WARN] Skipping PDF {path.name}: {exc}", file=sys.stderr)


# ── Page builders ─────────────────────────────────────────────────

def section_title_figure(title: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(A4_W, A4_H))
    ax.text(0.5, 0.5, title, ha="center", va="center",
            fontsize=30, fontweight="bold", transform=ax.transAxes)
    ax.axis("off")
    return fig


# A4 landscape in PDF points (1 pt = 1/72 inch)
_A4_W_PT = 841.89
_A4_H_PT = 595.28
_MARGIN_PT = 25.0  # ~0.35 inch


def _add_svg_vector(writer: pypdf.PdfWriter, path: Path) -> None:
    """Convert SVG → PDF (true vector, no rasterisation) and add as an A4 page."""
    pdf_bytes = cairosvg.svg2pdf(url=str(path))
    reader = pypdf.PdfReader(io.BytesIO(pdf_bytes))
    page = reader.pages[0]

    svg_x0 = float(page.mediabox.left)
    svg_y0 = float(page.mediabox.bottom)
    svg_w  = float(page.mediabox.width)
    svg_h  = float(page.mediabox.height)

    if svg_w == 0 or svg_h == 0:
        return

    avail_w = _A4_W_PT - 2 * _MARGIN_PT
    avail_h = _A4_H_PT - 2 * _MARGIN_PT
    scale   = min(avail_w / svg_w, avail_h / svg_h)

    # Centre on the page; shift by origin in case mediabox doesn't start at (0,0)
    tx = _MARGIN_PT + (avail_w - svg_w * scale) / 2 - scale * svg_x0
    ty = _MARGIN_PT + (avail_h - svg_h * scale) / 2 - scale * svg_y0

    page.add_transformation((scale, 0, 0, scale, tx, ty))
    page.mediabox.lower_left  = (0, 0)
    page.mediabox.upper_right = (_A4_W_PT, _A4_H_PT)

    writer.add_page(page)


def image_figure(path: Path) -> plt.Figure:
    """Render a PNG image on an A4 landscape matplotlib figure (raster fallback)."""
    img = Image.open(str(path))
    fig = plt.figure(figsize=(A4_W, A4_H))
    ax = fig.add_axes([0.02, 0.06, 0.96, 0.90])
    ax.imshow(img)
    ax.axis("off")
    fig.text(0.5, 0.01, path.name, ha="center", va="bottom",
             fontsize=7, color="#666666")
    return fig


def _table_figure(df: pd.DataFrame, title: str, page_label: str = "",
                  rows_per_page: int = 22) -> plt.Figure:
    """Render a DataFrame as a single A4 landscape figure with a styled table."""
    n_rows, n_cols = df.shape

    fig = plt.figure(figsize=(A4_W, A4_H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    # Title — fixed position
    y = 0.96
    ax.text(0.5, y, title, ha="center", va="top",
            fontsize=13, fontweight="bold", transform=ax.transAxes, color=HEADER_BG)
    y -= 0.06
    # Always render the page-label slot (text only when non-empty) so that
    # table_top is identical on every page regardless of how full the page is.
    if page_label:
        ax.text(0.5, y, page_label, ha="center", va="top",
                fontsize=9, color="#7f8c8d", transform=ax.transAxes)
    y -= 0.04   # reserved unconditionally

    table_top = y - 0.01

    # Format cell values
    def fmt(v):
        if pd.isna(v):
            return "—"
        if isinstance(v, float):
            return f"{v:.4g}"
        return str(v)

    cell_vals = [[fmt(v) for v in row] for row in df.itertuples(index=False)]

    tbl = ax.table(
        cellText=cell_vals,
        colLabels=df.columns.tolist(),
        cellLoc="left",
        loc="upper center",
        bbox=[0.01, 0.02, 0.98, table_top - 0.02],
    )
    tbl.auto_set_font_size(False)
    font_size = max(5, min(9, int(72 / n_cols)))
    tbl.set_fontsize(font_size)
    tbl.auto_set_column_width(col=list(range(n_cols)))

    for j in range(n_cols):
        tbl[0, j].set_facecolor(HEADER_BG)
        tbl[0, j].set_text_props(color=HEADER_FG, fontweight="bold")

    for i in range(1, n_rows + 1):
        for j in range(n_cols):
            if i % 2 == 0:
                tbl[i, j].set_facecolor(ALT_ROW)

    # Fix every cell to the same height so a partial last page looks identical
    # to a full page.  Height is expressed in axes-coordinate units.
    fixed_row_h = (table_top - 0.02) / (rows_per_page + 1)  # +1 for header
    for cell in tbl.get_celld().values():
        cell.set_height(fixed_row_h)

    return fig


def dataframe_figures(df: pd.DataFrame, title: str, rows_per_page: int = 22):
    """Yield one figure per page of a (possibly large) DataFrame."""
    total = len(df)
    n_pages = max(1, -(-total // rows_per_page))
    for i in range(n_pages):
        chunk = df.iloc[i * rows_per_page : (i + 1) * rows_per_page]
        label = f"Page {i + 1} of {n_pages}" if n_pages > 1 else ""
        yield _table_figure(chunk, title, label, rows_per_page)


def study_params_figure() -> plt.Figure:
    """Render rarefaction depths and (optionally) DADA2 truncation lengths."""
    fig, ax = plt.subplots(figsize=(A4_W, A4_H))
    ax.axis("off")

    ax.text(0.5, 0.93, "Study Parameters", ha="center", va="top",
            fontsize=22, fontweight="bold", transform=ax.transAxes, color=HEADER_BG)
    ax.plot([0.05, 0.95], [0.875, 0.875],
            transform=ax.transAxes, color=HEADER_BG, linewidth=1.5)

    y = 0.82

    def section_hdr(text):
        nonlocal y
        ax.text(0.08, y, text, ha="left", va="top",
                fontsize=14, fontweight="bold", transform=ax.transAxes, color=ACCENT)
        y -= 0.045
        ax.plot([0.08, 0.92], [y, y],
                transform=ax.transAxes, color=ACCENT, linewidth=0.8)
        y -= 0.025

    def item(label, value):
        nonlocal y
        ax.text(0.12, y, "•", ha="left", va="top",
                fontsize=12, transform=ax.transAxes, color="#95a5a6")
        ax.text(0.16, y, label, ha="left", va="top",
                fontsize=12, transform=ax.transAxes)
        ax.text(0.72, y, str(value), ha="left", va="top",
                fontsize=12, fontweight="bold", transform=ax.transAxes,
                color=VALUE_COLOR, fontfamily="monospace")
        y -= 0.058

    # Rarefaction depths
    section_hdr("Rarefaction Depths")
    with open(rarefy_base) as fh:
        depth = int(next(csv.DictReader(fh))["rarefaction_depth"])
    item("Non-phylogenetic  (Bray-Curtis / Jaccard)", f"{depth:,} reads")

    for db, path in zip(phylogenetic_reference_dbs, rarefy_dbs):
        with open(path) as fh:
            depth = int(next(csv.DictReader(fh))["rarefaction_depth"])
        item(f"{db}  (Faith PD / UniFrac)", f"{depth:,} reads")

    # Truncation lengths (NGS only)
    if trunc_len_csv and Path(trunc_len_csv).exists():
        y -= 0.03
        section_hdr("DADA2 Truncation Lengths")
        with open(trunc_len_csv) as fh:
            row = next(csv.DictReader(fh))
        item("Forward read (R1)", f"{row['trunc_len_f']} bp")
        item("Reverse read (R2)", f"{row['trunc_len_r']} bp")

    return fig


# Column rename map to tighten the permanova table display
_PERM_COL_RENAME = {
    "method name":          "method",
    "test statistic name":  "statistic name",
    "test statistic":       "statistic",
    "number of permutations": "permutations",
    "sample size":          "n",
    "p-value":              "p-value",
    "p_adjusted":           "p-adj",
}


def permanova_figures_for_db(path: Path, db_name: str):
    """Yield figures for one DB's PERMANOVA/PERMDISP TSV (global then pairwise)."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns=_PERM_COL_RENAME)

    for col in df.select_dtypes(include=["float64", "float32"]):
        df[col] = df[col].round(4)

    global_df = (
        df[df["scope"] == "global"]
        .drop(columns=["scope", "pair"], errors="ignore")
        .reset_index(drop=True)
    )
    pairwise_df = (
        df[df["scope"] == "pairwise"]
        .drop(columns=["scope"], errors="ignore")
        .reset_index(drop=True)
    )

    if not global_df.empty:
        yield from dataframe_figures(
            global_df, f"PERMANOVA & PERMDISP — {db_name}  (Global)"
        )
    if not pairwise_df.empty:
        yield from dataframe_figures(
            pairwise_df, f"PERMANOVA & PERMDISP — {db_name}  (Pairwise)"
        )


# ── Image section definitions ─────────────────────────────────────

_TAXA_LEVEL_RANK = {lvl: i for i, lvl in enumerate(
    ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
)}

def _sort_image_paths(section_title: str, paths) -> list:
    """Sort image paths; taxa barplots ordered Kingdom → Species, others alphabetically."""
    if section_title == "Taxonomic Composition":
        def _key(p):
            for lvl, rank in _TAXA_LEVEL_RANK.items():
                if f"taxa_barplot_{lvl}_" in p.name:
                    return (rank, p.name)
            return (len(_TAXA_LEVEL_RANK), p.name)
        return sorted(paths, key=_key)
    return sorted(paths)

SECTIONS = [
    ("Taxonomic Composition",           lambda p: "taxa_barplot"           in str(p)),
    ("Alpha Diversity",                  lambda p: "alpha_diversity"        in str(p)),
    ("Beta Diversity",                   lambda p: "beta_diversity"         in str(p)),
    ("Differential Abundance",           lambda p: "differential_abundance" in str(p)),
    ("Functional Prediction (PICRUSt2)", lambda p: "picrust2"              in str(p).lower()),
]

buckets = {title: [] for title, _ in SECTIONS}
for p in plot_inputs:
    for title, pred in SECTIONS:
        if pred(p):
            buckets[title].append(p)
            break

# ── Assemble PDF ──────────────────────────────────────────────────

writer = pypdf.PdfWriter()

# 1. Study Parameters
add_figure(writer, study_params_figure())

# 2. Metadata
meta_df = pd.read_csv(metadata_path, sep="\t")
for fig in dataframe_figures(meta_df, "Sample Metadata"):
    add_figure(writer, fig)

# 3. PERMANOVA / PERMDISP (global then pairwise, one set per DB)
for db, path in zip(reference_dbs, permanova_paths):
    for fig in permanova_figures_for_db(path, db):
        add_figure(writer, fig)

# 4. Image and PDF sections
for section_title, paths in buckets.items():
    image_paths = _sort_image_paths(section_title, [p for p in paths if p.suffix in (".png", ".svg")])
    pdf_paths   = sorted(p for p in paths if p.suffix == ".pdf")
    if not image_paths and not pdf_paths:
        continue
    add_figure(writer, section_title_figure(section_title))
    for p in image_paths:
        try:
            if p.suffix.lower() == ".svg":
                _add_svg_vector(writer, p)
            else:
                add_figure(writer, image_figure(p))
        except Exception as exc:
            print(f"[WARN] Skipping {p.name}: {exc}", file=sys.stderr)
    for p in pdf_paths:
        add_pdf_file(writer, p)

with open(out_pdf, "wb") as fh:
    writer.write(fh)

print(f"[REPORT] Visualisations PDF written → {out_pdf}", file=sys.stderr)
