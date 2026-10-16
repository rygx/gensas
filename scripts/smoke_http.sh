#!/usr/bin/env bash
# Checks a running GenSAS server that was started with a token: /version, authentication, a check, a run with its
# result, and the web client page (when the server serves it).
#
# Usage:  bash scripts/smoke_http.sh http://127.0.0.1:8765 TOKEN
set -euo pipefail

BASE="${1:?base url}"
TOKEN="${2:?token}"

for _ in $(seq 1 60); do
  if curl -fsS "${BASE}/version" >/dev/null 2>&1; then break; fi
  sleep 1
done
curl -fsS "${BASE}/version" | python3 -c "import sys, json; v = json.load(sys.stdin); assert v['authRequired'] and v['api'] == 1, v; print('version ok:', v['server'], v['engine']['version'])"

code=$(curl -s -o /dev/null -w '%{http_code}' "${BASE}/examples")
[ "${code}" = "401" ] || { echo "expected 401 without the token, got ${code}"; exit 1; }

AUTH=(-H "Authorization: Bearer ${TOKEN}" -H "Content-Type: application/json")
SRC='model D\n  parameter Real k = 1;\n  Real x(start = 1);\nequation\n  der(x) = -k*x;\nend D;\n'
curl -fsS "${AUTH[@]}" -d "{\"source\": \"${SRC}\"}" "${BASE}/check" | python3 -c "import sys, json; r = json.load(sys.stdin); assert r['ok'], r; print('check ok:', r['info']['name'])"

RUN=$(curl -fsS "${AUTH[@]}" -d "{\"model\": {\"source\": \"${SRC}\"}, \"solver\": {\"endTime\": 1}}" "${BASE}/runs" | python3 -c "import sys, json; print(json.load(sys.stdin)['id'])")
for _ in $(seq 1 60); do
  STATE=$(curl -fsS -H "Authorization: Bearer ${TOKEN}" "${BASE}/runs/${RUN}" | python3 -c "import sys, json; print(json.load(sys.stdin)['state'])")
  case "${STATE}" in completed) break ;; failed | cancelled) echo "run ${STATE}"; exit 1 ;; esac
  sleep 1
done
[ "${STATE}" = "completed" ] || { echo "run did not finish (${STATE})"; exit 1; }
curl -fsS -H "Authorization: Bearer ${TOKEN}" "${BASE}/runs/${RUN}/result" | python3 -c "
import sys, json, math
r = json.load(sys.stdin)
x = r['variables']['x'][-1]
assert abs(x - math.exp(-1)) < 1e-4, x
print('run ok: x(1) =', x)"

curl -fsS "${BASE}/" | grep -qi "<div id=\"app\"\|<title>" || { echo "the web client is not served"; exit 1; }
echo "client page ok"

