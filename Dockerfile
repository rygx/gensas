# GenSAS server with the web client, in one image.
#
#   docker build -t gensas .
#   docker run --rm -p 8764:8764 -e GENSAS_TOKEN=change-me gensas
#
# then open http://localhost:8764/ . See doc/Deployment.md for the LAN and Google Cloud set-ups.

# ---- the web client ---------------------------------------------------------------------
FROM node:22-bookworm-slim AS web
WORKDIR /src
COPY api ./api
WORKDIR /src/web
COPY web/package.json web/package-lock.json ./
RUN npm ci --no-audit --no-fund
COPY web/ ./
RUN npm run build

# ---- the engine (libgensas.so) -----------------------------------------------------------
FROM python:3.12-slim-bookworm AS engine
RUN apt-get update \
 && apt-get install -y --no-install-recommends g++ libarmadillo-dev libsuperlu-dev libjsoncpp-dev libblas-dev liblapack-dev \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY include ./include
COPY sas ./sas
COPY io ./io
COPY util ./util
COPY capi ./capi
COPY scripts/build_shared_lib.sh scripts/solver_sources.txt ./scripts/
# (the sed makes the build work from a Windows checkout whose files still have CRLF line endings)
RUN sed -i 's/\r$//' scripts/build_shared_lib.sh scripts/solver_sources.txt && bash scripts/build_shared_lib.sh /out

# ---- the image that runs ----------------------------------------------------------------
FROM python:3.12-slim-bookworm
RUN apt-get update \
 && apt-get install -y --no-install-recommends libjsoncpp25 libsuperlu5 libblas3 liblapack3 \
 && rm -rf /var/lib/apt/lists/* \
 && useradd --system --uid 10001 --no-create-home gensas
WORKDIR /app
COPY server/requirements.txt ./server/requirements.txt
RUN pip install --no-cache-dir -r server/requirements.txt
COPY python ./python
COPY server/gensas_server ./server/gensas_server
COPY --from=engine /out/libgensas.so ./lib/libgensas.so
COPY --from=web /src/web/dist ./web

ENV PYTHONPATH=/app/python:/app/server \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    GENSAS_LIBRARY=/app/lib/libgensas.so \
    GENSAS_STATIC_DIR=/app/web \
    GENSAS_HOST=0.0.0.0
# Results are written below the temporary directory, so the container also runs with --read-only --tmpfs /tmp.
USER gensas
EXPOSE 8764
# The server listens on $PORT (set by Cloud Run), or $GENSAS_PORT, or 8764.
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s CMD python -c "import os, urllib.request as u; u.urlopen('http://127.0.0.1:%s/version' % (os.environ.get('PORT') or os.environ.get('GENSAS_PORT') or '8764'), timeout=4)"
CMD ["python", "-m", "gensas_server"]
