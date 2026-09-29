"""Print a report to PDF with every panel view shown and labelled (the HTML is not modified).

The HTML panels show one view at a time (Species / Genus and above); the PDF
shows them all, each headed by its view, and drops the toggle buttons. Each barplot and cladogram
pair is stacked at full width, and only small units (an image, a caption, a table row, a short
finding) are kept from splitting across pages.

Usage, with the host python3 (Playwright, chromium headless shell):
    python3 report_to_pdf.py <report.html> <report.pdf> [<report.html> <report.pdf> ...]
"""
import sys
from playwright.sync_api import sync_playwright

CSS = """
.panel .view{display:block !important; margin-bottom:14px}
.seg, .seg-all{display:none !important}
.panel .view > .view-label{font:600 12px/1.4 var(--mono, monospace); color:#555;
  text-transform:uppercase; letter-spacing:.06em; margin:10px 0 4px; break-after:avoid}
/* Barplot above cladogram, each at full page width: side by side they are unreadable on A4. */
.pair{display:block !important}
/* Each panel image at most about half a page, so a barplot and its cladogram can share one;
   at full width each fills a page and cannot be split, leaving a gap before it. A tall image
   (a barplot with many rows) is the exception: capped at half a page it is too narrow to
   read, so it keeps the full width (class set in JS below). */
.pair img{width:auto !important; max-width:100% !important; max-height:128mm !important;
  margin:0 auto 10px !important}
.pair img.tall{width:100% !important; max-height:240mm !important}
/* A figure taller than a page cannot stay in one piece; clamp each image to a page. */
figure img{max-height:240mm; object-fit:contain}
/* Keep only small units together. Keeping a whole multi-view panel or a long table together
   pushes it to the next page and leaves the rest of the current one blank. */
img, figcaption, tr, .finding, .note, .tile{break-inside:avoid}
h2, h3, .view-label, thead, .panel-head, .rowhead{break-after:avoid}
/* The small-multiple grids (alpha diversity tiles) stay two to a row: the page is narrower than
   the report's 760px breakpoint, which would otherwise blow each tile up to a full page. */
.grid2{grid-template-columns:1fr 1fr !important}
/* On screen wide tables scroll sideways; on paper the overflow is clipped. Let text cells wrap
   (numbers and chips stay on one line), tighten padding, and print at SCALE so every table fits. */
.scroll{overflow:visible !important}
th, td{white-space:normal !important; padding-right:8px !important}
td.n, th.n, .chip{white-space:nowrap !important}
table{font-size:.8rem !important}
table.long td.q{min-width:0 !important}
/* A status chip with a qualifier ("q < 0.05, not robust") may wrap; the others are one word. */
table.long td:last-child .chip{white-space:normal !important; display:inline-block}
"""
SCALE = 0.85
CONTENT_MM = 186   # A4 width minus the 12 mm side margins
JS = """
() => {
  const names = {sp: 'Species', gn: 'Genus and above', pa: 'Per-arm fit', j: 'Joint fit'};
  document.querySelectorAll('.pair img').forEach(i => {
    if (i.naturalHeight > 0.8 * i.naturalWidth) i.classList.add('tall');
  });
  document.querySelectorAll('.panel .view').forEach(v => {
    const parts = [...v.classList].filter(c => c in names).map(c => names[c]);
    const order = ['Per-arm fit', 'Joint fit', 'Species', 'Genus and above'];
    parts.sort((a, b) => order.indexOf(a) - order.indexOf(b));
    const d = document.createElement('div'); d.className = 'view-label';
    d.textContent = parts.join(' · '); v.prepend(d);
  });
}
"""
with sync_playwright() as p:
    b = p.chromium.launch()
    for src, out in zip(sys.argv[1::2], sys.argv[2::2]):
        # Lay the page out at the width it will print at, so screen breakpoints behave as on paper.
        pg = b.new_page(viewport={"width": int(CONTENT_MM / 25.4 * 96 / SCALE), "height": 1200})
        pg.emulate_media(color_scheme="light")
        pg.goto("file://" + src); pg.wait_for_load_state("networkidle")
        pg.add_style_tag(content=CSS); pg.evaluate(JS)
        pg.emulate_media(media="print", color_scheme="light")
        broken = pg.evaluate("[...document.images].filter(i=>!i.complete||i.naturalWidth===0).length")
        clipped = pg.evaluate("[...document.querySelectorAll('.scroll')].filter(e => (e.firstElementChild||e).scrollWidth > e.clientWidth + 1).length")
        if broken or clipped:
            raise SystemExit(f"{src}: {broken} broken image(s), {clipped} table(s) wider than the page")
        n = pg.evaluate("document.querySelectorAll('.panel .view').length")
        pg.pdf(path=out, format="A4", print_background=True, scale=SCALE,
               margin={"top": "14mm", "bottom": "14mm", "left": "12mm", "right": "12mm"})
        print(f"{out}: {n} panel views shown, broken images {broken}")
        pg.close()
    b.close()
