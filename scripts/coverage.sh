#!/usr/bin/env bash
# Measures how much of the code the tests run: line and branch coverage of the C++ library
# (sas/, io/, util/, capi/) by the C++ unit tests and the Python tests together, and line coverage
# of the Python wrapper (python/gensas).
#
# It uses the same Bazel-free build as scripts/local_solver_test.sh, so it covers what that script
# builds: the Modelica path, results, C API and case files. The power-flow code and the PSAT
# reader are tested through Bazel only and are not part of the numbers.
#
# Usage (from the repository root):
#   bash scripts/coverage.sh
#
# Prerequisites: those of local_solver_test.sh, plus  pip install gcovr coverage
#
# Output (COVERAGE_DIR, default /tmp/gensas-coverage):
#   summary.md          totals and the per-file table (also written to the CI job summary)
#   cpp/index.html      browsable per-line report of the C++ code
#   cpp.json            gcovr JSON (per-line data)
#   python/index.html   browsable report of the Python wrapper
#
# Variables: those of local_solver_test.sh, and
#   COVERAGE_DIR       output directory
#   SHARED_LINK_LIBS   libraries for linking the shared library the Python tests load (default: LINK_LIBS)
#   SKIP_PYTHON=1      only measure the C++ code

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

OUT="${COVERAGE_DIR:-/tmp/gensas-coverage}"
rm -rf "${OUT}"
mkdir -p "${OUT}"

export BUILD_DIR="${OUT}/build"
# -fPIC: the same objects also go into the shared library for the Python tests.
export EXTRA_CXXFLAGS="--coverage -fprofile-update=atomic -O0 -fPIC ${EXTRA_CXXFLAGS:-}"
export EXTRA_LDFLAGS="--coverage"

echo "== C++ tests"
bash scripts/local_solver_test.sh

if [ "${SKIP_PYTHON:-0}" != "1" ]; then
	echo "== Python tests"
	SHARED_LINK_LIBS="${SHARED_LINK_LIBS:-${LINK_LIBS:--ljsoncpp -lsuperlu -llapack -lblas -lpthread -lm}}"
	# The objects built above carry coverage counters; the shared library is only the test subject here.
	${CXX:-g++} -shared --coverage -o "${OUT}/libgensas.so" "${BUILD_DIR}"/obj/*.o ${SHARED_LINK_LIBS}
	GENSAS_LIBRARY="${OUT}/libgensas.so" python3 -m coverage run --data-file="${OUT}/python.coverage" \
		--include='*/python/gensas/*' -m unittest discover -s python/tests
	python3 -m coverage report --data-file="${OUT}/python.coverage" | tee "${OUT}/python.txt"
	python3 -m coverage html --data-file="${OUT}/python.coverage" -d "${OUT}/python" >/dev/null
	python3 -m coverage json --data-file="${OUT}/python.coverage" -o "${OUT}/python.json" >/dev/null
fi

echo "== C++ coverage"
mkdir -p "${OUT}/cpp"
GCOVR_ARGS=(--root . --object-directory "${BUILD_DIR}" --filter '(sas|io|util|capi)/' --exclude '.*/test/.*'
	--exclude-throw-branches --exclude-unreachable-branches --gcov-ignore-parse-errors=negative_hits.warn_once_per_file)
gcovr "${GCOVR_ARGS[@]}" --json "${OUT}/cpp.json" --html-details "${OUT}/cpp/index.html" --txt "${OUT}/cpp.txt"
cat "${OUT}/cpp.txt"

python3 scripts/coverage_summary.py "${OUT}" > "${OUT}/summary.md"
echo
cat "${OUT}/summary.md"
echo
echo "Reports: ${OUT}/cpp/index.html and ${OUT}/python/index.html"
