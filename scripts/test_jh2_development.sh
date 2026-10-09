#!/usr/bin/env bash

set -euo pipefail

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
build_dir=${PERIDIGM_JH2_BUILD_DIR:-"${source_dir}-build-jh2"}
image=${PERIDIGM_JH2_IMAGE:-"localhost/peridigm-jh2-dev:f25c492"}

podman run --rm --userns=keep-id \
  --volume "${source_dir}:/workspace/source:ro" \
  --volume "${build_dir}:/workspace/build" \
  --workdir /workspace/build \
  "${image}" \
  ctest --output-on-failure -R 'utJH2MaterialPoint|utPeridigm_JH2CorrespondenceMaterial'
