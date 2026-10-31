#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image_name="${MUTSEQR_BIOC_DEVEL_IMAGE:-mutseqr-bioc-devel:local}"
dockerfile="${repo_root}/docker/bioc-devel.Dockerfile"

build_image() {
  docker build \
    --pull \
    --progress plain \
    --file "${dockerfile}" \
    --tag "${image_name}" \
    "${repo_root}"
}

ensure_image() {
  if ! docker image inspect "${image_name}" >/dev/null 2>&1; then
    echo "Bioconductor devel image '${image_name}' is not available; building it now."
    build_image
  fi
}

run_container() {
  docker run \
    --rm \
    --init \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "${repo_root}:/workspace" \
    --workdir /workspace \
    "${image_name}" \
    "$@"
}

run_info() {
  run_container Rscript -e 'cat("R ", R.version.string, "\nBioconductor ", as.character(BiocManager::version()), "\n", sep = "")'
}

run_test() {
  run_container Rscript scripts/run-tests.R
}

run_check() {
  run_container Rscript scripts/run-check.R
}

run_bioccheck() {
  run_container Rscript -e 'BiocCheck::BiocCheck(".", `quit-with-status` = TRUE)'
}

command="${1:-}"

case "${command}" in
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
  validate)
    ensure_image
    run_info
    run_test
    run_check
    run_bioccheck
    ;;
  shell)
    ensure_image
    docker run \
      --rm \
      --init \
      --interactive \
      --tty \
      --user "$(id -u):$(id -g)" \
      --env HOME=/tmp \
      --volume "${repo_root}:/workspace" \
      --workdir /workspace \
      "${image_name}" \
      bash
    ;;
  *)
    echo "Usage: $0 {build|info|test|check|bioccheck|validate|shell}" >&2
    exit 2
    ;;
esac
