#!/usr/bin/env bash
# Builds the Docker image (unless IMAGE is given), starts it and checks that the server, the engine and the web client
# work: /version, a check, a run with its result, the client page, and the limits of the container.
#
# Usage (from the repository root):  bash scripts/docker_smoke_test.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
IMAGE="${IMAGE:-gensas-smoke}"
NAME="gensas-smoke-$$"
PORT_HOST="${PORT_HOST:-18764}"
TOKEN="smoke-test-token"

if [ -z "${SKIP_BUILD:-}" ]; then
  docker build -t "${IMAGE}" .
fi

cleanup() { docker rm -f "${NAME}" >/dev/null 2>&1 || true; }
trap cleanup EXIT

# read-only root file system, as the deployment guide recommends
docker run -d --name "${NAME}" --read-only --tmpfs /tmp -p "${PORT_HOST}:8764" -e GENSAS_TOKEN="${TOKEN}" -e GENSAS_RATE_LIMIT=60 "${IMAGE}" >/dev/null
BASE="http://127.0.0.1:${PORT_HOST}"

bash scripts/smoke_http.sh "${BASE}" "${TOKEN}"

[ "$(docker exec "${NAME}" id -u)" != "0" ] || { echo "the container runs as root"; exit 1; }
echo "runs as a non-root user"
echo "docker smoke test passed"
