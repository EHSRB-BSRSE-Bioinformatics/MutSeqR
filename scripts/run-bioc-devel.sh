#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image_name="${MUTSEQR_BIOC_DEVEL_IMAGE:-mutseqr-bioc-devel:local}"
dockerfile="${repo_root}/docker/bioc-devel.Dockerfile"

# Portable on Linux and macOS; detects stale images after dependency/tool edits.
build_inputs="$(cksum < "$dockerfile")/$(cksum < "${repo_root}/DESCRIPTION")"

build_image() {
  docker build \
    --pull \
    --progress plain \
    --file "${dockerfile}" \
    --tag "${image_name}" \
    --label "org.mutseqr.validation-inputs=${build_inputs}" \
    "${repo_root}"
}

ensure_image() {
  if ! docker info >/dev/null 2>&1; then
    echo "Docker is unavailable. Start Docker Engine/Desktop and retry." >&2
    exit 1
  fi
  local installed_inputs
  installed_inputs=$(docker image inspect --format '{{index .Config.Labels "org.mutseqr.validation-inputs"}}' "$image_name" 2>/dev/null || true)
  if [[ "$installed_inputs" != "$build_inputs" ]]; then
    echo "Building Bioconductor validation image (missing or dependency inputs changed)."
    build_image
  fi
}

run_container() {
  mkdir -p "${repo_root}/.validation"
  mkdir -p "${repo_root}/.cache/ExperimentHub"
  docker run \
    --rm \
    --init \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --env EXPERIMENT_HUB_CACHE=/cache/ExperimentHub \
    --env MUTSEQR_RESULTS_DIR=/results \
    --env MUTSEQR_CHECK_DIR=/results/check \
    --volume "${repo_root}:/source:ro" \
    --volume "${repo_root}/.cache/ExperimentHub:/cache/ExperimentHub" \
    --volume "${repo_root}/.validation:/results" \
    --entrypoint /bin/bash \
    "${image_name}" \
    /source/scripts/container-entrypoint.sh "$@"
}

run_info() {
  docker image inspect --format 'Validation image: {{.Id}}' "$image_name"
  run_container Rscript -e 'cat(R.version.string, "\nBioconductor ", as.character(BiocManager::version()), "\n", sep = "")'
}

run_test() {
  run_container Rscript scripts/run-tests.R
}

run_check() {
  run_container Rscript scripts/run-check.R
}

run_bioccheck() {
  run_container Rscript scripts/run-bioccheck.R
}

command="${1:-}"

case "${command}" in
  fingerprint)
    echo "$build_inputs"
    ;;
  build)
    build_image
    ;;
  info)
    ensure_image
    run_info
    ;;
  test)
    ensure_image
    run_test
    ;;
  check)
    ensure_image
    run_check
    ;;
  bioccheck)
    ensure_image
    run_bioccheck
    ;;
  validate|release-check)
    ensure_image
    mkdir -p "${repo_root}/.validation"
    {
      run_info
      run_container Rscript scripts/run-validation.R "$command"
    } 2>&1 | tee "${repo_root}/.validation/validation.log"
    ;;
  shell)
    ensure_image
    # Interactive shells use the same isolated workspace as validation.
    mkdir -p "${repo_root}/.validation"
    docker run --rm --init --interactive --tty \
      --user "$(id -u):$(id -g)" --env HOME=/tmp \
      --env MUTSEQR_RESULTS_DIR=/results --env MUTSEQR_CHECK_DIR=/results/check \
      --volume "${repo_root}:/source:ro" \
      --volume "${repo_root}/.validation:/results" \
      --entrypoint /bin/bash "$image_name" \
      /source/scripts/container-entrypoint.sh bash
    ;;
  *)
    echo "Usage: $0 {build|info|test|check|bioccheck|validate|release-check|shell|fingerprint}" >&2
    exit 2
    ;;
esac
