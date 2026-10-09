#!/usr/bin/env bash

set -euo pipefail

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
build_dir=${PERIDIGM_JH2_BUILD_DIR:-"${source_dir}-build-jh2"}
image=${PERIDIGM_JH2_IMAGE:-"localhost/peridigm-jh2-dev:f25c492"}

mkdir -p "${build_dir}"

podman run --rm --userns=keep-id \
  --volume "${source_dir}:/workspace/source:ro" \
  --volume "${build_dir}:/workspace/build" \
  --workdir /workspace/build \
  "${image}" \
  cmake \
  -S /workspace/source \
  -B /workspace/build \
  -D CMAKE_BUILD_TYPE:STRING=Debug \
  -D CMAKE_CXX_STANDARD:STRING=17 \
  -D CMAKE_CXX_STANDARD_REQUIRED:BOOL=ON \
  -D CMAKE_C_COMPILER:FILEPATH=/usr/bin/mpicc \
  -D CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/mpicxx \
  -D CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/mpifort \
  -D TRILINOS_DIR:PATH=/usr/local/trilinos

podman run --rm --userns=keep-id \
  --volume "${source_dir}:/workspace/source:ro" \
  --volume "${build_dir}:/workspace/build" \
  --workdir /workspace/build \
  "${image}" \
  cmake --build /workspace/build --parallel 4
