#!/usr/bin/env bash
# Builds and runs the solver unit tests (sas/test/*_test.cpp) without Bazel.
#
# Why: a Bazel build compiles OpenBLAS, SuperLU, ... from source, which takes many minutes. The
# solver tests only need Armadillo, SuperLU, jsoncpp and a BLAS/LAPACK, so with those installed on
# the system a change-compile-test loop takes seconds. CI (Bazel) stays the authority.
#
# Usage (from the repository root):
#   bash scripts/local_solver_test.sh                 # all tests
#   bash scripts/local_solver_test.sh hybrid_sim_test # selected tests
#
# Prerequisites, Debian/Ubuntu:
#   sudo apt install g++ libarmadillo-dev libsuperlu-dev libjsoncpp-dev libblas-dev liblapack-dev
# Other systems: install the same libraries and point the variables below at them.
#
# Variables (defaults in parentheses):
#   BUILD_DIR          object and binary directory (/tmp/gensas-local)
#   CXX                compiler (g++)
#   ARMADILLO_INCLUDE  directory containing armadillo (/usr/include)
#   SUPERLU_INCLUDE    directory containing slu_ddefs.h (/usr/include/superlu)
#   JSONCPP_INCLUDE    directory containing json/json.h (/usr/include/jsoncpp)
#   LINK_LIBS          libraries to link (-ljsoncpp -lsuperlu -llapack -lblas -lpthread -lm)
#   JOBS               parallel compile jobs (number of CPUs)
#   EXTRA_CXXFLAGS     added to every compile (scripts/coverage.sh passes --coverage -O0)
#   EXTRA_LDFLAGS      added to every link
#

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

BUILD_DIR="${BUILD_DIR:-/tmp/gensas-local}"
CXX="${CXX:-g++}"
ARMADILLO_INCLUDE="${ARMADILLO_INCLUDE:-/usr/include}"
SUPERLU_INCLUDE="${SUPERLU_INCLUDE:-/usr/include/superlu}"
JSONCPP_INCLUDE="${JSONCPP_INCLUDE:-/usr/include/jsoncpp}"
LINK_LIBS="${LINK_LIBS:--ljsoncpp -lsuperlu -llapack -lblas -lpthread -lm}"
JOBS="${JOBS:-$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null) || echo 4)}"

CXXFLAGS="-std=c++14 -O1 -w -DARMA_DONT_USE_WRAPPER -I. \
-I${ARMADILLO_INCLUDE} -I${SUPERLU_INCLUDE} -I${JSONCPP_INCLUDE} ${EXTRA_CXXFLAGS:-}"

# Library sources the solver tests depend on (no power flow, no PSAT/MAT input); the list is shared with
# scripts/build_shared_lib.sh.
SOURCES="$(tr '\n' ' ' < scripts/solver_sources.txt)"

mkdir -p "${BUILD_DIR}/obj" "${BUILD_DIR}/bin"

# Compile the library objects that are missing or older than their source or any header.
newest_header=$(find sas util io -name '*.h' -o -name '*.txt' | xargs ls -t | head -n 1)
compile() {
  local src="$1"
  local obj="${BUILD_DIR}/obj/$(basename "${src}").o"
  if [ ! -f "${obj}" ] || [ "${src}.cpp" -nt "${obj}" ] || [ "${newest_header}" -nt "${obj}" ]; then
    echo "CXX ${src}.cpp"
    ${CXX} ${CXXFLAGS} -c "${src}.cpp" -o "${obj}"
  fi
}
export -f compile
export BUILD_DIR CXX CXXFLAGS newest_header
echo ${SOURCES} | tr ' ' '\n' | xargs -P "${JOBS}" -I{} bash -c 'compile {}'

# Tests live in sas/test and capi/test (.cpp, or .c for the C header check).
test_source() {
  for ext in cpp c; do
    for dir in sas/test capi/test; do
      if [ -f "${dir}/$1.${ext}" ]; then echo "${dir}/$1.${ext}"; return; fi
    done
  done
}
if [ "$#" -gt 0 ]; then
  TESTS="$*"
else
  TESTS="$(ls sas/test/*_test.cpp capi/test/*_test.c* | sed 's|.*/||; s/\.[a-z]*$//')"
fi

FAILED=0
for t in ${TESTS}; do
  echo "=== ${t} ==="
  src="$(test_source "${t}")"
  case "${src}" in
    *.c)
      ${CC:-gcc} -std=c99 -Wall -I. -c "${src}" -o "${BUILD_DIR}/${t}.o"
      ${CXX} "${BUILD_DIR}/${t}.o" "${BUILD_DIR}"/obj/*.o ${LINK_LIBS} ${EXTRA_LDFLAGS:-} -o "${BUILD_DIR}/bin/${t}"
      ;;
    *)
      ${CXX} ${CXXFLAGS} "${src}" "${BUILD_DIR}"/obj/*.o ${LINK_LIBS} ${EXTRA_LDFLAGS:-} -o "${BUILD_DIR}/bin/${t}"
      ;;
  esac
  # Run in the build directory: some tests write result files into the working directory.
  if ! (cd "${BUILD_DIR}" && "./bin/${t}" > "${t}.log" 2>&1); then
    echo "FAILED: ${t} (see ${BUILD_DIR}/${t}.log)"
    FAILED=1
  else
    grep -E 'PASSED' "${BUILD_DIR}/${t}.log" | tail -n 1
  fi
done

exit ${FAILED}
