#!/usr/bin/env bash
set -euo pipefail

# Check an isolated copy, including uncommitted edits, without rewriting sources
# or leaving generated plots in the host worktree.
workspace=$(mktemp -d /tmp/mutseqr-workspace-XXXXXX)
cleanup() {
  # Some BiocCheck versions ignore checkDir and write beside the source tree.
  # Preserve those reports on success and failure before removing the workspace.
  if [[ -d "$workspace/MutSeqR.BiocCheck" && -n "${MUTSEQR_RESULTS_DIR:-}" ]]; then
    cp -R "$workspace/MutSeqR.BiocCheck" "$MUTSEQR_RESULTS_DIR/"
  fi
  rm -rf "$workspace"
}
trap cleanup EXIT
# BiocCheck requires the source directory to match DESCRIPTION Package.
mkdir "$workspace/MutSeqR"
tar --directory /source --exclude='./.git' --exclude='./.pixi' \
    --exclude='./.cache' \
    --exclude='./.validation' -cf - . | tar --directory "$workspace/MutSeqR" -xf -
cd "$workspace/MutSeqR"
"$@"
