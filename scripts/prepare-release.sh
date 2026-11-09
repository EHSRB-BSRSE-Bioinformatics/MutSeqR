#!/usr/bin/env bash
# Prepare the next devel version, then validate that exact source tree.
# Review and commit the resulting metadata yourself; this script never pushes.
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ $# -ne 0 ]]; then
  echo "Usage: pixi run -e container prepare-release (no automatic commit or push)" >&2
  exit 2
fi
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Working tree is dirty; commit or stash your changes first." >&2
  exit 1
fi
old_ver=$(awk '/^Version:/ {print $2; exit}' DESCRIPTION)
if [[ ! "$old_ver" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "DESCRIPTION Version must have the form x.y.z" >&2
  exit 1
fi
z=${old_ver##*.}
new_ver="${old_ver%.*}.$((10#$z + 1))"
news_top=$(awk '/^# / {print; exit}' NEWS.md)
if [[ "$news_top" != '# MutSeqR (development version)' && "$news_top" != "# MutSeqR $new_ver" ]]; then
  echo "Start NEWS.md with '# MutSeqR (development version)' or '# MutSeqR $new_ver'." >&2
  exit 1
fi

# POSIX awk and temporary files also work on macOS (unlike GNU sed -i).
# Stage and check both files before changing either tracked file.
release_tmp=$(mktemp -d)
trap 'rm -rf "$release_tmp"' EXIT
awk -v version="$new_ver" '/^Version:/ {$0 = "Version: " version} {print}' DESCRIPTION > "$release_tmp/DESCRIPTION"
awk -v version="$new_ver" '$0 == "# MutSeqR (development version)" && !changed { $0 = "# MutSeqR " version; changed = 1 } {print}' NEWS.md > "$release_tmp/NEWS.md"
prepared_ver=$(awk '/^Version:/ {print $2; exit}' "$release_tmp/DESCRIPTION")
prepared_news_top=$(awk '/^# / {print; exit}' "$release_tmp/NEWS.md")
if [[ "$prepared_ver" != "$new_ver" || "$prepared_news_top" != "# MutSeqR $new_ver" ]]; then
  echo "Prepared DESCRIPTION and NEWS.md must both match version $new_ver." >&2
  exit 1
fi
mv "$release_tmp/DESCRIPTION" DESCRIPTION
mv "$release_tmp/NEWS.md" NEWS.md

echo "Prepared $old_ver -> $new_ver; validating the final version."
bash scripts/run-bioc-devel.sh release-check
echo "Review DESCRIPTION and NEWS.md, then commit and push when ready."
