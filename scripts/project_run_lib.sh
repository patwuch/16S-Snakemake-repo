# Helpers for the numbered project analyses under main/<study>/additional_comparison_scripts/run_*/.
# Sourced by each analysis's env.sh, never run on its own. No pipeline rule uses this file.
#
# What it provides:
#   in_sif NAME CMD...        run CMD in containers/NAME.sif the way the pipeline does
#                             (apptainer exec --home <repo>, from the repo root)
#   light LABEL CMD...        run CMD (a quick step: always reruns)
#   heavy LABEL CMD...        run CMD unless every path in the OUTS array exists and is newer
#                             than every path in the INS array; FORCE=1 runs it regardless
#   check_metadata META STUDY_META
#                             stop unless STUDY_META still has the md5 recorded in META.source
#   record_metadata META STUDY_META
#                             archive META under metadata_versions/ and record STUDY_META's md5

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FORCE="${FORCE:-0}"
DRY_RUN="${DRY_RUN:-0}"

die() { echo "ERROR: $*" >&2; exit 1; }

step() { printf '\n== %s\n' "$*"; }

in_sif() {
    local sif="$REPO/containers/$1.sif"; shift
    [ -e "$sif" ] || die "missing container $sif"
    (cd "$REPO" && apptainer exec --home "$REPO" "$sif" "$@")
}

light() {
    local label="$1"; shift
    if [ "$DRY_RUN" = 1 ]; then
        echo "   would run $label"
        return 0
    fi
    echo "   run $label"
    "$@"
}

heavy() {
    local label="$1"; shift
    if [ "$FORCE" != 1 ] && _fresh; then
        echo "   skip $label (outputs newer than inputs)"
        return 0
    fi
    if [ "$DRY_RUN" = 1 ]; then
        echo "   would run $label"
        return 0
    fi
    echo "   run $label"
    "$@"
}

_fresh() {
    local o i
    for o in "${OUTS[@]}"; do
        [ -e "$o" ] || return 1
        for i in "${INS[@]}"; do
            [ -e "$i" ] || die "missing input $i"
            [ "$i" -nt "$o" ] && return 1
        done
    done
    return 0
}

check_metadata() {
    local meta="$1" study="$2" src="$1.source"
    [ -e "$meta" ] || die "missing $meta. Write it by hand from $study, then: run_all.sh --record-metadata"
    [ -e "$src" ] || die "no $src. Check $meta against $study, then: run_all.sh --record-metadata"
    local want have
    want="$(awk '{print $1; exit}' "$src")"
    have="$(md5sum "$study" | cut -d' ' -f1)"
    [ "$want" = "$have" ] || die "$study changed since $(basename "$meta") was made ($(cut -d' ' -f3- "$src")).
Review the change. If the analysis metadata needs it, edit $meta.
Either way, then run: run_all.sh --record-metadata"
}

record_metadata() {
    local meta="$1" study="$2"
    [ -e "$meta" ] || die "missing $meta"
    local dir name last
    dir="$(dirname "$meta")/metadata_versions"
    name="$(basename "$meta" .tsv)"
    mkdir -p "$dir"
    last="$(ls -1 "$dir"/"$name".*.tsv 2>/dev/null | tail -n1 || true)"
    if [ -z "$last" ] || ! cmp -s "$meta" "$last"; then
        cp -p "$meta" "$dir/$name.$(date +%Y-%m-%dT%H%M%S).tsv"
        echo "archived $(basename "$meta") under $dir"
    fi
    echo "$(md5sum "$study" | cut -d' ' -f1)  $(basename "$study")  recorded $(date +%Y-%m-%d)" > "$meta.source"
    echo "recorded $study for $(basename "$meta")"
}

# run_steps DIR "$@": the shared run_all.sh body. Options:
#   --force             rerun the heavy steps even when their outputs are newer than their inputs
#   --dry-run           report what would run; run nothing
#   --from NN           start at step NN
#   --only NN           run step NN alone
#   --record-metadata   record the current study metadata for this analysis's metadata, then stop
run_steps() {
    local dir="$1"; shift
    local from=0 only=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --force) export FORCE=1 ;;
            --dry-run) export DRY_RUN=1 ;;
            --from) from="$((10#$2))"; shift ;;
            --only) only="$((10#$2))"; shift ;;
            --record-metadata) export RECORD_METADATA=1 ;;
            *) die "unknown option $1" ;;
        esac
        shift
    done
    if [ "${RECORD_METADATA:-0}" = 1 ]; then
        bash -c "source '$dir/env.sh' --no-check && record_metadata \"\$META\" \"\$STUDY_META\""
        return
    fi
    local s n
    for s in "$dir"/[0-9][0-9]_*.sh; do
        n="$((10#$(basename "$s" | cut -c1-2)))"
        [ -n "$only" ] && [ "$n" != "$only" ] && continue
        [ "$n" -lt "$from" ] && continue
        step "$(basename "$s" .sh)"
        bash "$s"
    done
}
