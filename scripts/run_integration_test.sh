#!/usr/bin/env bash
# Integration test: run power flow on d_70k_070.mat and verify:
#   1. The solver converges to alpha=1
#   2. "Computation time:" is printed (clean exit after computation)
#   3. Output .mat file is created and non-empty
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
echo "Input : resources/psat_mat/d_70k_070.mat"
echo "Output: ${OUTPUT_MAT}"
echo ""

# Build first (exits non-zero on build failure, aborting the test).
bazel build //app:app

echo ""

# Run the binary directly — bypasses the Bazel sandbox so the output
# .mat file can be written to an arbitrary filesystem path.
# Pass a relative output name; the app prepends CWD, yielding REPO_DIR/OUTPUT_NAME.
set +e
APP_OUTPUT=$("${APP_BIN}" \
  -p \
  -f "${INPUT_DIR}/resources/psat_mat/d_70k_070.mat" \
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

# ── Check 4: output .mat file exists and is non-empty ───────────────────────
if [ -f "${OUTPUT_MAT}" ]; then
  MAT_SIZE=$(wc -c < "${OUTPUT_MAT}")
  if [ "${MAT_SIZE}" -gt 0 ]; then
    echo "PASS: output .mat created (${MAT_SIZE} bytes)"
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
rm -f "${REPO_DIR}/${HYBRID_JSON}"
set +e
MO_OUTPUT=$("${APP_BIN}" \
  -g -m file \
  -i "${INPUT_DIR}/resources/mofile/test_bouncing_ball.mo" \
  -t 1.2 \
  -j "${HYBRID_JSON}" \
  2>&1)
MO_EXIT=$?
set -e

echo "=== GenSAS Modelica hybrid integration test ==="
echo "${MO_OUTPUT}" | grep -E "^Event at|Json output|rror" || true
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
rm -f "${REPO_DIR}/${HYBRID_JSON}"

# ── Summary ─────────────────────────────────────────────────────────────────
echo ""
if [ "${PASS}" -eq 1 ]; then
  echo "=== All checks passed ==="
  exit 0
else
  echo "=== One or more checks FAILED ==="
  exit 1
fi
