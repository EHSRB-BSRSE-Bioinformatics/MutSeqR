#!/usr/bin/env bash
# Exercise release preparation in disposable repositories without R or Docker.
set -euo pipefail
cd "$(dirname "$0")/.."
source_root=$PWD
test_root=$(mktemp -d)
trap 'rm -rf "$test_root"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

fixture() {
  local name=$1 heading=$2 validation_status=$3
  repo="$test_root/$name"
  mkdir -p "$repo/scripts"
  cp "$source_root/scripts/prepare-release.sh" "$repo/scripts/"
  printf 'Package: MutSeqR\nVersion: 1.1.9\n' > "$repo/DESCRIPTION"
  printf '%s\n\n- Release note.\n\n# MutSeqR 1.1.8\n' "$heading" > "$repo/NEWS.md"
  cat > "$repo/scripts/run-bioc-devel.sh" <<EOF
set -euo pipefail
[[ \$# -eq 1 && \$1 == release-check ]]
grep -qx 'Version: 1.1.10' DESCRIPTION
[[ \$(awk '/^# / {print; exit}' NEWS.md) == '# MutSeqR 1.1.10' ]]
touch validation-ran
exit $validation_status
EOF
  git -C "$repo" init -q
  git -C "$repo" add .
  git -C "$repo" -c user.name='Release test' -c user.email='release@example.invalid' \
    -c core.hooksPath=/dev/null commit -qm fixture
}

run_prepare() {
  status=0
  bash "$repo/scripts/prepare-release.sh" > "$test_root/output" 2>&1 || status=$?
}

for heading in '# MutSeqR (development version)' '# MutSeqR 1.1.10'; do
  fixture valid "$heading" 0
  run_prepare
  [[ $status -eq 0 ]] || { cat "$test_root/output" >&2; fail "valid heading rejected"; }
  [[ -f "$repo/validation-ran" ]] || fail "prepared metadata was not validated"
  grep -qx 'Version: 1.1.10' "$repo/DESCRIPTION" || fail "patch version not incremented"
  [[ $(awk '/^# / {print; exit}' "$repo/NEWS.md") == '# MutSeqR 1.1.10' ]] || fail "NEWS not finalized"
  grep -qx '# MutSeqR 1.1.8' "$repo/NEWS.md" || fail "historical heading changed"
  rm -rf "$repo"
done

for heading in '# MutSeqR 1.1.9' '# MutSeqR (development version) typo' 'No heading'; do
  fixture invalid "$heading" 0
  run_prepare
  [[ $status -eq 1 ]] || fail "invalid heading accepted"
  [[ -z $(git -C "$repo" status --porcelain) ]] || fail "invalid heading modified the worktree"
  [[ ! -e "$repo/validation-ran" ]] || fail "validation ran for invalid input"
  rm -rf "$repo"
done

fixture validation-failure '# MutSeqR (development version)' 7
run_prepare
[[ $status -eq 7 ]] || fail "release-check exit status was not preserved"
[[ -f "$repo/validation-ran" ]] || fail "release-check did not run"
grep -qx 'Version: 1.1.10' "$repo/DESCRIPTION" || fail "prepared version not retained for review"
[[ $(awk '/^# / {print; exit}' "$repo/NEWS.md") == '# MutSeqR 1.1.10' ]] || fail "prepared NEWS not retained for review"

echo 'Release-preparation regression tests passed.'