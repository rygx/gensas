#!/usr/bin/env bash
# Measures how much of the code the tests run: line and branch coverage of the C++ library
# (sas/, io/, util/, capi/) by the C++ unit tests and the Python tests together, and line coverage
# of the Python wrapper (python/gensas), of the server (server/gensas_server) and of the web client (web/src).
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
#   SKIP_PYTHON=1      only measure the C++ code (also skips the server)
#   SKIP_SERVER=1      do not measure the server (it needs `pip install -r server/requirements-dev.txt`)
#   SKIP_WEB=1         do not measure the web client (it needs `npm ci` in web/)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

OUT="${COVERAGE_DIR:-/tmp/gensas-coverage}"
rm -rf "${OUT}"
mkdir -p "${OUT}"

export BUILD_DIR="${OUT}/build"
# -fPIC: the same objects also go into the shared library for the Python tests.
# Counters are updated atomically so that the multi-threaded tests do not corrupt them (GCC 11 and later;
# Clang does not know the option).
ATOMIC=""
if ${CXX:-g++} --version 2>/dev/null | grep -qi 'gcc\|g++' && ! ${CXX:-g++} --version 2>/dev/null | grep -qi clang; then
	ATOMIC="-fprofile-update=atomic"
fi
export EXTRA_CXXFLAGS="--coverage ${ATOMIC} -O0 -fPIC ${EXTRA_CXXFLAGS:-}"
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

	# The server runs checks and simulations in worker processes; their code (worker.py) is exercised by the
	# tests but not measured here.
	if [ "${SKIP_SERVER:-0}" != "1" ] && python3 -c "import fastapi, httpx, yaml, jsonschema, openapi_spec_validator, pytest" 2>/dev/null; then
		echo "== Server tests"
		(cd server && GENSAS_LIBRARY="${OUT}/libgensas.so" python3 -m coverage run --data-file="${OUT}/server.coverage" \
			--include='*/server/gensas_server/*' --omit='*/worker.py' -m pytest tests -q)
		python3 -m coverage report --data-file="${OUT}/server.coverage" | tee "${OUT}/server.txt"
		python3 -m coverage html --data-file="${OUT}/server.coverage" -d "${OUT}/server" >/dev/null
		python3 -m coverage json --data-file="${OUT}/server.coverage" -o "${OUT}/server.json" >/dev/null
	fi
fi

if [ "${SKIP_WEB:-0}" != "1" ] && [ -d web/node_modules ]; then
	echo "== Web client tests"
	(cd web && npx vitest run --coverage --coverage.reportsDirectory="${OUT}/web")
fi

echo "== C++ coverage"
mkdir -p "${OUT}/cpp"
GCOVR_ARGS=(--root . --object-directory "${BUILD_DIR}" --filter '(sas|io|util|capi)/' --exclude '.*/test/.*'
	--exclude-throw-branches --exclude-unreachable-branches --gcov-ignore-parse-errors=negative_hits.warn_once_per_file)
python3 -m gcovr "${GCOVR_ARGS[@]}" --json "${OUT}/cpp.json" --html-details "${OUT}/cpp/index.html" --txt "${OUT}/cpp.txt"
cat "${OUT}/cpp.txt"

python3 scripts/coverage_summary.py "${OUT}" > "${OUT}/summary.md"
echo
cat "${OUT}/summary.md"
echo
echo "Reports: ${OUT}/cpp/index.html, ${OUT}/python/index.html, ${OUT}/server/index.html, ${OUT}/web/index.html"
