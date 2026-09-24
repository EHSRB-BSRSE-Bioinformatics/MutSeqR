#!/usr/bin/env bash
# MutSeqR release preparation (Bioconductor "devel" push workflow).
#
# Steps:
#   1. regenerate docs (roxygen)
#   2. bump the z version in DESCRIPTION (in-place, format-preserving)
#   3. check NEWS.md is current (warns/stops if the top entry is stale)
#   4. run unit tests
#   5. commit docs + DESCRIPTION (you review NEWS.md and the rest)
#
# After this script finishes you still do, by hand:
#   git push origin main
#   git push upstream main:devel     # upstream = git.bioconductor.org
#
# Usage: bash scripts/prepare-release.sh [--push]
#   --push  also pushes to origin main (not to Bioconductor)

set -euo pipefail
cd "$(dirname "$0")/.."

do_push=0
[[ "${1:-}" == "--push" ]] && do_push=1

if [[ -n "$(git status --porcelain)" ]]; then
  echo "Working tree is dirty; commit or stash your changes first:"
  git status --short
  exit 1
fi

## 1. documentation ---------------------------------------------------
echo "==> 1/5 Regenerating documentation (roxygen)"
Rscript -e 'devtools::document()'

## 2. version bump ----------------------------------------------------
echo "==> 2/5 Bumping z version in DESCRIPTION"
old_ver="$(grep -m1 '^Version:' DESCRIPTION | awk '{print $2}')"
z="${old_ver##*.}"
new_ver="${old_ver%.*}.$(( z + 1 ))"
sed -i "s/^Version: ${old_ver}\$/Version: ${new_ver}/" DESCRIPTION
echo "    ${old_ver} -> ${new_ver}"

## 3. NEWS.md check ---------------------------------------------------
echo "==> 3/5 Checking NEWS.md"
NEWS_TOP="$(grep -m1 '^# ' NEWS.md || true)"
case "${NEWS_TOP}" in
  *" ${new_ver}"*)
    echo "    OK: NEWS.md has an entry for ${new_ver}" ;;
  *" ${old_ver}"*)
    echo "    REMINDER: top NEWS.md entry is still ${old_ver}."
    echo "    Consider adding a '# MutSeqR ${new_ver}' entry." ;;
  *)
    echo "    ERROR: top NEWS.md entry ('${NEWS_TOP}') matches neither"
    echo "    the old (${old_ver}) nor the new (${new_ver}) version."
    echo "    Update NEWS.md before pushing."
    exit 1 ;;
esac

## 4. tests ------------------------------------------------------------
echo "==> 4/5 Running unit tests"
Rscript -e 'devtools::load_all(quiet = TRUE); devtools::test()'

## 5. commit ------------------------------------------------------------
echo "==> 5/5 Committing version bump + regenerated docs"
git add DESCRIPTION NEWS.md NAMESPACE man/
git commit -m "chore: version bump ${new_ver}" --allow-empty

if [[ "${do_push}" -eq 1 ]]; then
  git push origin HEAD
  echo "Pushed to origin. Next (needs Bioc SSH key): git push upstream main:devel"
else
  echo
  echo "Done. Review the commit, then:"
  echo "  git push origin main"
  echo "  git push upstream main:devel   # upstream = git.bioconductor.org:packages/MutSeqR.git"
fi
