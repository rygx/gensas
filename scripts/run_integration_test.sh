#!/usr/bin/env bash
# Integration test: convert the PSAT file d_70k_070.mat to a case file, run power flow on it and verify:
#   1. The solver converges to alpha=1
#   2. "Computation time:" is printed (clean exit after computation)
#   3. Output .mat file is created and non-empty
#   and that the converged state matches reference values (also for the small d_014_ind_zip1 case)
# and run the Modelica bouncing-ball example (-g mode) and verify:
#   4. It exits cleanly, reports its 3 state events and writes a JSON result with events
#
# Usage: bash scripts/run_integration_test.sh
#   Must be run from the repository root.
#
# Note: the binary is invoked directly (not via "bazel run") to avoid
# the Bazel darwin sandbox restricting write access to arbitrary paths.

set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# The app prepends GetCurrentWorkingDir() to the -o argument, so we pass
# a bare filename and let the app place it in the CWD (REPO_DIR).
OUTPUT_NAME="gensas_integration_test_output.mat"
OUTPUT_MAT="${REPO_DIR}/${OUTPUT_NAME}"

cd "${REPO_DIR}"
rm -f "${OUTPUT_MAT}"

# Windows (MSYS2 bash): the binary is app.exe, and the native executable needs a
# Windows-style path (C:/...) rather than an MSYS-style one (/c/...).
if command -v cygpath >/dev/null 2>&1; then
  INPUT_DIR="$(cygpath -m "${REPO_DIR}")"
else
  INPUT_DIR="${REPO_DIR}"
fi
APP_BIN="${REPO_DIR}/bazel-bin/app/app"
if [ -f "${APP_BIN}.exe" ]; then
  APP_BIN="${APP_BIN}.exe"
fi

echo "=== GenSAS power flow integration test ==="
echo "Input : resources/psat_mat/d_70k_070.mat (converted to a case file)"
echo "Output: ${OUTPUT_MAT}"
echo ""

# Build first (exits non-zero on build failure, aborting the test).
bazel build //app:app

echo ""

# The solver reads GenSAS case files; convert the PSAT data first (python3 with numpy, scipy and h5py).
PY="$(command -v python3 || command -v python || true)"
CASE_DIR="gensas_integration_cases"
rm -rf "${REPO_DIR}/${CASE_DIR}"
mkdir -p "${REPO_DIR}/${CASE_DIR}"
for name in d_70k_070 d_014_ind_zip1; do
  if [ -z "${PY}" ] || ! "${PY}" scripts/psat_to_case.py "resources/psat_mat/${name}.mat" "${CASE_DIR}/${name}.case" >/dev/null; then
    echo "FAIL: ${name}: could not convert the PSAT file to a case file (python3 with numpy, scipy and h5py required)"
    exit 1
  fi
done

# Run the binary directly — bypasses the Bazel sandbox so the output
# .mat file can be written to an arbitrary filesystem path.
# Pass a relative output name; the app prepends CWD, yielding REPO_DIR/OUTPUT_NAME.
set +e
APP_OUTPUT=$("${APP_BIN}" \
  -p \
  -f "${INPUT_DIR}/${CASE_DIR}/d_70k_070.case" \
  -s 0.5 -l 28 -d 1e-5 \
  -o "${OUTPUT_NAME}" \
  2>&1)
APP_EXIT=$?
set -e

echo "${APP_OUTPUT}"
echo ""

PASS=1

# ── Check 1: process exited cleanly ────────────────────────────────────────
if [ "${APP_EXIT}" -ne 0 ]; then
  echo "FAIL: app exited with code ${APP_EXIT} (crash or error)"
  PASS=0
else
  echo "PASS: clean exit (exit code 0)"
fi

# ── Check 2: solver converged to alpha=1 ────────────────────────────────────
if echo "${APP_OUTPUT}" | grep -q "Step=1,"; then
  echo "PASS: solver converged to alpha=1"
else
  echo "FAIL: 'Step=1,' not found — solver did not reach alpha=1"
  PASS=0
fi

# ── Check 3: timing line present (confirms no crash after computation) ───────
if echo "${APP_OUTPUT}" | grep -q "Computation time:"; then
  echo "PASS: 'Computation time:' found"
else
  echo "FAIL: 'Computation time:' not found — likely crashed before clean exit"
  PASS=0
fi

# Compares the "Final state" fingerprint printed by the app with reference values.
# The converged state is unique, so it agrees across platforms and BLAS builds to far better than
# the 1e-6 relative tolerance used here (the step sequence does not).
check_fingerprint() {
  local label="$1" output="$2" want_n="$3" want_norm="$4" want_sum="$5"
  local line
  line=$(echo "${output}" | grep "^Final state:" | tail -n 1)
  if [ -z "${line}" ]; then
    echo "FAIL: ${label}: no 'Final state' line"
    PASS=0
    return
  fi
  local n norm sum
  n=$(echo "${line}" | sed -E 's/.*n=([0-9]+),.*/\1/')
  norm=$(echo "${line}" | sed -E 's/.*norm=([-0-9.eE+]+),.*/\1/')
  sum=$(echo "${line}" | sed -E 's/.*sum=([-0-9.eE+]+).*/\1/')
  if [ "${n}" = "${want_n}" ] && awk -v a="${norm}" -v b="${want_norm}" -v c="${sum}" -v d="${want_sum}" \
    'function rel(x, y) { return (x - y) / (y < 0 ? -y : y) } BEGIN { r1 = rel(a, b); r2 = rel(c, d); if (r1 < 0) r1 = -r1; if (r2 < 0) r2 = -r2; exit !(r1 < 1e-6 && r2 < 1e-6) }'; then
    echo "PASS: ${label}: final state matches reference (n=${n}, norm=${norm}, sum=${sum})"
  else
    echo "FAIL: ${label}: final state differs; got n=${n} norm=${norm} sum=${sum}, want n=${want_n} norm=${want_norm} sum=${want_sum}"
    PASS=0
  fi
}

# ── Check 3b: converged state matches the reference ──────────────────────────
check_fingerprint "d_70k_070" "${APP_OUTPUT}" 490000 282.7520234 69169.68943

# ── Small power-flow case ────────────────────────────────────────────────────
set +e
SMALL_OUTPUT=$("${APP_BIN}" \
  -p \
  -f "${INPUT_DIR}/${CASE_DIR}/d_014_ind_zip1.case" \
  -s 0.5 -l 28 -d 1e-5 \
  -o "gensas_integration_test_small.mat" \
  2>&1)
SMALL_EXIT=$?
set -e
rm -f "${REPO_DIR}/gensas_integration_test_small.mat"
if [ "${SMALL_EXIT}" -ne 0 ]; then
  echo "FAIL: small power-flow case exited with code ${SMALL_EXIT}"
  PASS=0
else
  check_fingerprint "d_014_ind_zip1" "${SMALL_OUTPUT}" 153 4.061145417 11.64662517
fi

# SciPy (an independent implementation) must be able to read the MAT files we write.
check_mat_with_scipy() {
  local label="$1" file="$2" code="$3"
  if [ -n "${PY}" ] && "${PY}" -c "import sys; from scipy.io import loadmat; m = loadmat(sys.argv[1]); ${code}" "${file}"; then
    echo "PASS: ${label}: SciPy reads the MAT file"
  else
    echo "FAIL: ${label}: SciPy could not read the MAT file or its content is wrong"
    PASS=0
  fi
}

# ── Check 4: output .mat file exists and is non-empty ───────────────────────
if [ -f "${OUTPUT_MAT}" ]; then
  MAT_SIZE=$(wc -c < "${OUTPUT_MAT}")
  if [ "${MAT_SIZE}" -gt 0 ]; then
    echo "PASS: output .mat created (${MAT_SIZE} bytes)"
    check_mat_with_scipy "power-flow result" "${INPUT_DIR}/${OUTPUT_NAME}" "assert m['s'].shape == (490000, 1)"
  else
    echo "FAIL: output .mat file is empty"
    PASS=0
  fi
else
  echo "FAIL: output .mat file was not created at ${OUTPUT_MAT}"
  PASS=0
fi

# ── Modelica (hybrid) command line ──────────────────────────────────────────
# Bouncing ball: three state events (reinit of the velocity) before t = 1.2.
HYBRID_JSON="gensas_integration_test_hybrid.json"
HYBRID_CSV="gensas_integration_test_hybrid.csv"
HYBRID_EVENTS_CSV="gensas_integration_test_hybrid.events.csv"
HYBRID_MAT="gensas_integration_test_hybrid.mat"
rm -f "${REPO_DIR}/${HYBRID_JSON}" "${REPO_DIR}/${HYBRID_CSV}" "${REPO_DIR}/${HYBRID_EVENTS_CSV}" "${REPO_DIR}/${HYBRID_MAT}"
set +e
MO_OUTPUT=$("${APP_BIN}" \
  -g -m file \
  -i "${INPUT_DIR}/resources/mofile/test_bouncing_ball.mo" \
  -t 1.2 \
  -j "${HYBRID_JSON}" \
  -c "${HYBRID_CSV}" \
  -o "${HYBRID_MAT}" \
  2>&1)
MO_EXIT=$?
set -e

echo "=== GenSAS Modelica hybrid integration test ==="
echo "${MO_OUTPUT}" | grep -E "^Event at|output written|rror" || true
echo ""

if [ "${MO_EXIT}" -ne 0 ]; then
  echo "FAIL: Modelica run exited with code ${MO_EXIT}"
  PASS=0
else
  echo "PASS: Modelica run exited cleanly"
fi

MO_EVENTS=$(echo "${MO_OUTPUT}" | grep -c "^Event at" || true)
if [ "${MO_EVENTS}" -eq 3 ]; then
  echo "PASS: 3 bounce events reported"
else
  echo "FAIL: expected 3 bounce events, found ${MO_EVENTS}"
  PASS=0
fi

if [ -s "${REPO_DIR}/${HYBRID_JSON}" ] && grep -q '"events"' "${REPO_DIR}/${HYBRID_JSON}"; then
  echo "PASS: JSON result with events written"
else
  echo "FAIL: JSON result with events was not written to ${HYBRID_JSON}"
  PASS=0
fi
if [ -s "${REPO_DIR}/${HYBRID_CSV}" ] && head -n 1 "${REPO_DIR}/${HYBRID_CSV}" | grep -q '^time,' && [ -s "${REPO_DIR}/${HYBRID_EVENTS_CSV}" ]; then
  echo "PASS: CSV result and events written"
else
  echo "FAIL: CSV result or events file was not written (${HYBRID_CSV})"
  PASS=0
fi

# A MAT v5 file starts with the text header "MATLAB 5.0 MAT-file".
if [ -s "${REPO_DIR}/${HYBRID_MAT}" ] && head -c 19 "${REPO_DIR}/${HYBRID_MAT}" | grep -q "MATLAB 5.0 MAT-file"; then
  echo "PASS: MAT v5 result written"
  check_mat_with_scipy "Modelica result" "${INPUT_DIR}/${HYBRID_MAT}" "assert m['s'].shape[0] == len(m['names']) and m['s'].shape[1] == m['t'].shape[1] and m['event_t'].size == 3"
else
  echo "FAIL: MAT v5 result was not written to ${HYBRID_MAT}"
  PASS=0
fi
rm -f "${REPO_DIR}/${HYBRID_JSON}" "${REPO_DIR}/${HYBRID_CSV}" "${REPO_DIR}/${HYBRID_EVENTS_CSV}" "${REPO_DIR}/${HYBRID_MAT}"

# ── Modelica job file (gensas run) ──────────────────────────────────────────
# Paths in the job file are relative to the job file's directory (the repository root here).
JOB_FILE="gensas_integration_job.json"
JOB_JSON="gensas_integration_job_out.json"
cat > "${REPO_DIR}/${JOB_FILE}" <<'JOB'
{
  "gensas": {"format": "job", "version": 1},
  "kind": "modelica",
  "model": {"file": "resources/mofile/test_bouncing_ball.mo"},
  "solver": {"endTime": 1.2},
  "output": {"json": "gensas_integration_job_out.json"}
}
JOB
rm -f "${REPO_DIR}/${JOB_JSON}"
set +e
JOB_OUTPUT=$("${APP_BIN}" run "${INPUT_DIR}/${JOB_FILE}" 2>&1)
JOB_EXIT=$?
BAD_OUTPUT=$("${APP_BIN}" run "${INPUT_DIR}/${JOB_FILE}" --set solver.endTim=3 2>&1)
BAD_EXIT=$?
set -e
echo "${JOB_OUTPUT}" | grep -E "^Event at|output written|rror" || true
if [ "${JOB_EXIT}" -eq 0 ] && [ "$(echo "${JOB_OUTPUT}" | grep -c '^Event at' || true)" -eq 3 ] && [ -s "${REPO_DIR}/${JOB_JSON}" ]; then
  echo "PASS: gensas run produced the bounce events and the JSON result"
else
  echo "FAIL: gensas run (exit ${JOB_EXIT})"
  PASS=0
fi
if [ "${BAD_EXIT}" -ne 0 ] && echo "${BAD_OUTPUT}" | grep -q "unknown key solver.endTim"; then
  echo "PASS: gensas run rejects an unknown key"
else
  echo "FAIL: gensas run accepted an unknown key (exit ${BAD_EXIT})"
  PASS=0
fi
rm -f "${REPO_DIR}/${JOB_FILE}" "${REPO_DIR}/${JOB_JSON}"

rm -rf "${REPO_DIR}/${CASE_DIR}"

# ── Summary ─────────────────────────────────────────────────────────────────
echo ""
if [ "${PASS}" -eq 1 ]; then
  echo "=== All checks passed ==="
  exit 0
else
  echo "=== One or more checks FAILED ==="
  exit 1
fi
