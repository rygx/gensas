# Deployment

The [server](Server.md) and the [web client](Web.md) are packaged in one Docker image: the server answers the API and also serves the client page. This page shows how to run it on your own computer, for other computers at home, and on Google Cloud.

## The image

```bash
docker build -t gensas .
docker run --rm -p 8765:8765 -e GENSAS_TOKEN=change-me gensas
```

Open `http://localhost:8765/`, press **Server** and enter the token. The image runs as an unprivileged user and needs no writable file system except `/tmp`, so this is the stricter way to run it:

```bash
docker run -d --name gensas --restart unless-stopped --read-only --tmpfs /tmp \
  --memory 1g --cpus 1 --pids-limit 256 --cap-drop ALL --security-opt no-new-privileges \
  -p 8765:8765 -e GENSAS_TOKEN=change-me gensas
```

`bash scripts/docker_smoke_test.sh` builds the image and checks it (the Docker workflow does the same on every push). The settings are the environment variables of the [server](Server.md#settings); `PORT` is read too, because Cloud Run sets it.

## On your computer and at home

* Only you, one machine, no Docker: follow [Server](Server.md) (`python -m gensas_server`), it listens on `127.0.0.1` only. The client from `web/` (`npm run dev`) or the image above both work.
* Other computers on your home network: run the image (or the server with `--host 0.0.0.0`) on one machine and open `http://<that machine's address>:8765/` from the others. **Always set a token** (`GENSAS_TOKEN`) when the port is reachable by other computers; everything except `/version` then needs it. Do not forward the port on your router: the connection has no encryption, so the token and the models travel in the clear. For access from outside, use a VPN (WireGuard, Tailscale) to your home network, or one of the public set-ups below, which have TLS.
* Each run is a separate process. Set `GENSAS_MAX_CONCURRENT_RUNS` to the number of CPU cores you want to give it (the default is 2).

## A public app on Google Cloud

A public server accepts models from anyone, so it needs limits more than a private one does. What the image already does: every check and run is a separate process with a time limit (`GENSAS_RUN_TIMEOUT`, 120 s), the size of the model text is limited, the number of runs waiting is limited (429), results are removed after an hour, and the process is unprivileged. Add:

* `GENSAS_RATE_LIMIT=30`: at most 30 checks and runs a minute for each client address (429 with `Retry-After` beyond that). Behind a proxy or load balancer also set `GENSAS_TRUST_PROXY=1` so that the address is read from `X-Forwarded-For` (only do this when the server cannot be reached except through your proxy, or the header can be forged).
* Small memory limits for the container, and `GENSAS_MAX_CONCURRENT_RUNS` no larger than the CPUs; fewer kept results (`GENSAS_MAX_KEPT_RUNS=20`, `GENSAS_RETENTION=900`).
* A **budget alert** in the Google Cloud console (Billing, Budgets & alerts), and a hard cap on the number of instances.
* No token for a public demo (a token in a page everyone can open is no secret), or a token you hand out to the people who may use it.

### Option 1: Cloud Run (nothing to maintain, pay for use)

Cloud Run runs the container behind HTTPS and a public address. Two properties of the server matter: it keeps runs and results in memory and on the local disk of one instance, so it must run as **one instance** (`--max-instances 1`); and a run goes on after the request that started it has returned, so the CPU must stay allocated (`--no-cpu-throttling`), otherwise a run stalls between polls.

```bash
gcloud config set project YOUR_PROJECT
gcloud services enable run.googleapis.com cloudbuild.googleapis.com artifactregistry.googleapis.com

# build the image with Cloud Build and deploy it (run from the repository root)
gcloud run deploy gensas --source . --region us-central1 --allow-unauthenticated \
  --max-instances 1 --concurrency 40 --memory 1Gi --cpu 1 --no-cpu-throttling --timeout 900 \
  --set-env-vars GENSAS_RATE_LIMIT=30,GENSAS_TRUST_PROXY=1,GENSAS_MAX_CONCURRENT_RUNS=1,GENSAS_MAX_KEPT_RUNS=20,GENSAS_RETENTION=900
```

The command prints the address of the service. The file system of Cloud Run lives in memory, so stored results count against `--memory`. With `--no-cpu-throttling` an instance is billed for the whole time it exists; add `--min-instances 0` (the default) to let it stop when nobody uses it, at the price of a few seconds of start-up for the first visitor. To require a token, add `GENSAS_TOKEN` (best from a Secret Manager secret: `--set-secrets GENSAS_TOKEN=NAME:latest`) and give it to the users.

### Option 2: one virtual machine (predictable cost, a little upkeep)

A small Compute Engine VM running the same image costs a fixed amount each month; the free `e2-micro` (1 GB memory, shared CPU, in `us-west1`, `us-central1` or `us-east1`) can run light use but has little network traffic included and is slow for long simulations, so `e2-small` or `e2-medium` is more comfortable. Install Docker on the VM, run the image as above with `--restart unless-stopped`, and put a TLS proxy in front of it; with [Caddy](https://caddyserver.com) and a domain name pointing at the VM the whole proxy configuration is

```
gensas.example.com {
    reverse_proxy 127.0.0.1:8765
}
```

(then publish the container on `127.0.0.1:8765` only, and set `GENSAS_TRUST_PROXY=1`). Open ports 80 and 443 in the VM's firewall rules and nothing else.

### Which one

For a hobby app that is used now and then, Cloud Run is the cheapest and needs no maintenance. For steady use or long simulations a VM costs the same every month whatever happens, and the work stays on one machine. Prices and free amounts change; check [Cloud Run pricing](https://cloud.google.com/run/pricing) and [Compute Engine pricing](https://cloud.google.com/compute/all-pricing) before deciding, and keep the budget alert on either way.

### Serving the client elsewhere

The client is a static site (`web/`, `npm run build`, the `dist/` folder) that talks to the API at the address typed under **Server**. It can be hosted on any static host (Firebase Hosting, Cloud Storage, GitHub Pages) with the server somewhere else; the server then needs `GENSAS_CORS_ORIGINS=https://your-client-host`. Serving both from the image avoids that.
