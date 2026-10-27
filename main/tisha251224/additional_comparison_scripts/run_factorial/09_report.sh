# Report checks, then the figure zip, standalone HTML, share zip and PDF.
# report.html is hand-edited and only read. Host python3 (Playwright for the PDF).
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "report outputs" python3 "$SCRIPTS/build_report_outputs.py"
