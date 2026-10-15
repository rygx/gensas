#!/usr/bin/env bash
# Builds and runs the solver unit tests (sas/test/*_test.cpp) without Bazel.
#
# Why: a Bazel build compiles OpenBLAS, HDF5, matio, ... from source, which takes many minutes. The
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
#
# matio is replaced by a stub (scripts/local_solver/matio_stub) for the library sources that include it;
# the solver and result tests do not need the real library. The MAT writer test is skipped.

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

CXXFLAGS="-std=c++14 -O1 -w -DARMA_DONT_USE_WRAPPER -I. -Iscripts/local_solver/matio_stub \
-I${ARMADILLO_INCLUDE} -I${SUPERLU_INCLUDE} -I${JSONCPP_INCLUDE}"

# Library sources the solver tests depend on (no power flow, no PSAT/MAT input).
SOURCES="sas/SasModel sas/SasConfig sas/SasInput sas/SasLexico sas/TokenScanner sas/ExprParser \
sas/ModelParser sas/ModelChecker sas/HybridModel sas/HybridSimulator sas/SasComputation sas/SasResult \
util/AbstractCheCalculator util/CheCompUtil util/CheState io/CheDataFormat io/CheIoUtil"

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
${CXX} ${CXXFLAGS} -c scripts/local_solver/matio_stub/matio_stub.cpp -o "${BUILD_DIR}/obj/matio_stub.o"

if [ "$#" -gt 0 ]; then
  TESTS="$*"
else
  # *_mat_test.cpp read MAT files back with matio and need the real library (Bazel builds them).
  TESTS="$(cd sas/test && ls *_test.cpp | grep -v '_mat_test' | sed 's/\.cpp$//')"
fi

FAILED=0
for t in ${TESTS}; do
  echo "=== ${t} ==="
  ${CXX} ${CXXFLAGS} "sas/test/${t}.cpp" "${BUILD_DIR}"/obj/*.o ${LINK_LIBS} -o "${BUILD_DIR}/bin/${t}"
  # Run in the build directory: some tests write result files into the working directory.
  if ! (cd "${BUILD_DIR}" && "./bin/${t}" > "${t}.log" 2>&1); then
    echo "FAILED: ${t} (see ${BUILD_DIR}/${t}.log)"
    FAILED=1
  else
    grep -E 'PASSED' "${BUILD_DIR}/${t}.log" | tail -n 1
  fi
done

exit ${FAILED}
