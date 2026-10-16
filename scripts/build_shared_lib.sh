#!/usr/bin/env bash
# Builds the shared library of the C API (libgensas.so, or .dylib on macOS) without Bazel, for the
# Python interface and the server tests. Same prerequisites and variables as scripts/local_solver_test.sh.
#
# Usage (from the repository root):
#   bash scripts/build_shared_lib.sh [OUTPUT_DIR]      # default /tmp/gensas-lib
#   GENSAS_LIBRARY=/tmp/gensas-lib/libgensas.so python3 -m unittest discover -s python/tests
#
# Variables: CXX, ARMADILLO_INCLUDE, SUPERLU_INCLUDE, JSONCPP_INCLUDE, LINK_LIBS, JOBS (see local_solver_test.sh).
#   SHARED_LINK_LIBS   libraries for the link (default: LINK_LIBS)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

OUT="${1:-/tmp/gensas-lib}"
CXX="${CXX:-g++}"
ARMADILLO_INCLUDE="${ARMADILLO_INCLUDE:-/usr/include}"
SUPERLU_INCLUDE="${SUPERLU_INCLUDE:-/usr/include/superlu}"
JSONCPP_INCLUDE="${JSONCPP_INCLUDE:-/usr/include/jsoncpp}"
LINK_LIBS="${LINK_LIBS:--ljsoncpp -lsuperlu -llapack -lblas -lpthread -lm}"
SHARED_LINK_LIBS="${SHARED_LINK_LIBS:-${LINK_LIBS}}"
JOBS="${JOBS:-$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null) || echo 4)}"
EXT="so"
if [ "$(uname)" = "Darwin" ]; then EXT="dylib"; fi

CXXFLAGS="-std=c++17 -O1 -w -fPIC -DARMA_DONT_USE_WRAPPER -I. -I${ARMADILLO_INCLUDE} -I${SUPERLU_INCLUDE} -I${JSONCPP_INCLUDE}"
mkdir -p "${OUT}/obj"

compile() {
  local src="$1"
  local obj="${OUT}/obj/$(basename "${src}").o"
  if [ ! -f "${obj}" ] || [ "${src}.cpp" -nt "${obj}" ]; then
    echo "CXX ${src}.cpp"
    ${CXX} ${CXXFLAGS} -c "${src}.cpp" -o "${obj}"
  fi
}
export -f compile
export OUT CXX CXXFLAGS
xargs -P "${JOBS}" -I{} bash -c 'compile {}' < scripts/solver_sources.txt

${CXX} -shared -o "${OUT}/libgensas.${EXT}" "${OUT}"/obj/*.o ${SHARED_LINK_LIBS}
echo "Built ${OUT}/libgensas.${EXT}"
