#!/usr/bin/env bash
# Report containers used by the workflow whose .sif no longer matches its
# sources (the .def and the files its %files section copies in). Exits 1 if
# any are stale; rebuilding is left to containers/build.sh on purpose, since
# unpinned rebuilds can silently change package versions.
#
# Fast path: a source is only inspected if it is newer than the .sif. Newer
# mtimes alone are not proof (git checkouts touch files), so those sources are
# compared against the copies baked into the image: the embedded deffile, and
# each %files destination path.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

stale=()
for sif in $(grep -rhoE 'containers/[A-Za-z0-9._-]+\.sif' Snakefile rules/ | sort -u); do
    def="${sif%.sif}.def"
    [ -f "$def" ] || continue
    if [ ! -f "$sif" ]; then
        stale+=("$(basename "$def" .def)  (no .sif built)")
        continue
    fi

    # "src dst" pairs from %files; the .def itself is compared via inspect.
    reasons=()
    if [ "$def" -nt "$sif" ] && ! cmp -s "$def" <(apptainer inspect --deffile "$sif" 2>/dev/null); then
        reasons+=("$def")
    fi
    while read -r src dst; do
        if [ "$src" -nt "$sif" ] && ! cmp -s "$src" <(apptainer exec "$sif" cat "$dst" 2>/dev/null); then
            reasons+=("$src")
        fi
    done < <(awk '/^%files/{f=1; next} /^%/{f=0} f && NF && $1 !~ /^#/ {print $1, $2}' "$def")

    if [ ${#reasons[@]} -gt 0 ]; then
        stale+=("$(basename "$def" .def)  (changed: ${reasons[*]})")
    fi
done

if [ ${#stale[@]} -gt 0 ]; then
    echo "Stale containers (sources differ from what is inside the .sif):" >&2
    printf '  %s\n' "${stale[@]}" >&2
    echo "Rebuild with: containers/build.sh --force NAME" >&2
    echo "Then rerun affected rules with -R; Snakemake will not do it on its own." >&2
    exit 1
fi
