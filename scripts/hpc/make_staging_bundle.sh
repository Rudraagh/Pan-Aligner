#!/bin/bash
# Build the HPC staging bundle on the laptop (Git Bash or WSL). Nothing is uploaded.
#
# Bundle = `git archive` of the committed HEAD (outputs/ excluded by default) plus the
# pristine upstream sources of the gitignored PanAligner/, minigraph/ and zlib/ checkouts
# at their recorded commits (no laptop .o files or Ubuntu-built binaries).
#
# With --wheelhouse DIR (from make_wheelhouse.sh), DIR is bundled as wheelhouse/ for the
# offline pbs_setup_venv.sh job and covered by the same checksums.
#
# Usage: bash scripts/hpc/make_staging_bundle.sh [--with-outputs] [--wheelhouse DIR] [OUT_DIR]
# Output (default hpc_staging/, gitignored):
#   panaligner_<commit>.tar.gz          the bundle (contains STAGED_FILES.sha256)
#   panaligner_<commit>.tar.gz.sha256   checksum to verify after upload
#   bundle_manifest.json                commit, tool commits, sizes, exclusions
#   pbs_env_probe.sh, pbs_unpack_bundle.sh  standalone jobs to upload beside the bundle
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WITH_OUTPUTS=0
WHEELHOUSE=""
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --with-outputs) WITH_OUTPUTS=1; shift ;;
    --wheelhouse) WHEELHOUSE="$(cd "${2:?--wheelhouse needs a directory}" && pwd)"; shift 2 ;;
    *) echo "ERROR: unknown option $1" >&2; exit 1 ;;
  esac
done
OUT_DIR="${1:-$ROOT_DIR/hpc_staging}"
TOOLS=(PanAligner minigraph zlib)
# Git links in this repository, so `git archive` alone leaves them as empty folders.
THIRD_PARTY=(tools/third_party/pbsim2 tools/third_party/minimap2)
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"  # POSIX form; GNU tar reads "C:/..." as a remote host

cd "$ROOT_DIR"
COMMIT="$(git rev-parse HEAD)"
NAME="panaligner_$(git rev-parse --short=12 HEAD)"
BUNDLE="$OUT_DIR/$NAME.tar.gz"
if [[ -e "$BUNDLE" ]]; then
  echo "ERROR: $BUNDLE already exists; remove it yourself or choose another OUT_DIR." >&2
  exit 1
fi

# The bundle is the commit, not the working tree: warn about uncommitted code.
uncommitted="$(git status --porcelain --untracked-files=no -- . ':(exclude)outputs' ':(exclude)data/queries')"
if [[ -n "$uncommitted" ]]; then
  echo "WARNING: these tracked changes are NOT in the bundle (commit them first if needed):" >&2
  echo "$uncommitted" >&2
fi

for tool in "${TOOLS[@]}"; do
  if ! git -C "$ROOT_DIR/$tool" rev-parse HEAD >/dev/null 2>&1; then
    echo "ERROR: $tool/ is not a git checkout; cannot export a pristine source tree." >&2
    exit 1
  fi
done
# Local edits to tool sources would silently be lost by exporting HEAD. zlib's Makefile and
# zconf.h are rewritten by its ./configure, so those two are expected and ignored.
for tool in PanAligner minigraph "${THIRD_PARTY[@]}"; do
  if ! git -C "$ROOT_DIR/$tool" rev-parse HEAD >/dev/null 2>&1; then
    echo "ERROR: $tool/ is not a git checkout; cannot export a pristine source tree." >&2
    exit 1
  fi
  if [[ -n "$(git -C "$ROOT_DIR/$tool" status --porcelain --untracked-files=no)" ]]; then
    echo "ERROR: $tool/ has modified tracked sources; refusing to stage HEAD without them." >&2
    git -C "$ROOT_DIR/$tool" status --short --untracked-files=no >&2
    exit 1
  fi
done

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK:?}"' EXIT

pathspec=(.)
[[ "$WITH_OUTPUTS" == "1" ]] || pathspec+=(':(exclude)outputs')
git archive --format=tar --prefix="$NAME/" HEAD -- "${pathspec[@]}" | tar -x -C "$WORK"

: > "$WORK/$NAME/STAGED_TOOLS.txt"
for tool in "${TOOLS[@]}" "${THIRD_PARTY[@]}"; do
  git -C "$ROOT_DIR/$tool" archive --format=tar --prefix="$NAME/$tool/" HEAD | tar -x -C "$WORK"
  printf '%s %s %s\n' "$tool" "$(git -C "$ROOT_DIR/$tool" rev-parse HEAD)" "$(git -C "$ROOT_DIR/$tool" remote get-url origin 2>/dev/null || echo unknown)" >> "$WORK/$NAME/STAGED_TOOLS.txt"
done

wheel_tags=""
if [[ -n "$WHEELHOUSE" ]]; then
  cp -r "$WHEELHOUSE" "$WORK/$NAME/wheelhouse"
  wheel_tags="$(cd "$WORK/$NAME/wheelhouse" && ls -d py* 2>/dev/null | tr '\n' ' ')"
fi

(cd "$WORK/$NAME" && find . -type f ! -name STAGED_FILES.sha256 -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > STAGED_FILES.sha256)
file_count="$(wc -l < "$WORK/$NAME/STAGED_FILES.sha256")"

(cd "$OUT_DIR" && tar -czf "$NAME.tar.gz" -C "$WORK" "$NAME" && sha256sum "$NAME.tar.gz" > "$NAME.tar.gz.sha256")
bundle_sha="$(cut -d' ' -f1 < "$BUNDLE.sha256")"
bundle_bytes="$(wc -c < "$BUNDLE" | tr -d ' ')"

# Committed copies, so the standalone jobs match the bundle exactly (and keep LF endings).
for job in pbs_env_probe.sh pbs_unpack_bundle.sh; do
  if ! git show "HEAD:scripts/hpc/$job" > "$OUT_DIR/$job" 2>/dev/null; then
    echo "WARNING: scripts/hpc/$job is not committed; copying the working-tree version." >&2
    tr -d '\r' < "$ROOT_DIR/scripts/hpc/$job" > "$OUT_DIR/$job"
  fi
done

tool_json=""
while read -r tool sha remote; do
  tool_json+="${tool_json:+, }\"$tool\": {\"commit\": \"$sha\", \"remote\": \"$remote\"}"
done < "$WORK/$NAME/STAGED_TOOLS.txt"
cat > "$OUT_DIR/bundle_manifest.json" <<EOF
{
  "bundle": "$NAME.tar.gz",
  "sha256": "$bundle_sha",
  "bytes": $bundle_bytes,
  "commit": "$COMMIT",
  "created": "$(date -Is)",
  "file_count": $file_count,
  "outputs_included": $([[ "$WITH_OUTPUTS" == "1" ]] && echo true || echo false),
  "uncommitted_changes_excluded": $([[ -n "$uncommitted" ]] && echo true || echo false),
  "wheelhouse": "${wheel_tags% }",
  "tools": {$tool_json}
}
EOF

cat "$OUT_DIR/bundle_manifest.json"
echo "Upload to the cluster: $NAME.tar.gz, $NAME.tar.gz.sha256, pbs_unpack_bundle.sh (and pbs_env_probe.sh)"
echo "from $OUT_DIR"
