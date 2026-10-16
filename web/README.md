# GenSAS web client

A single-page application (Preact, CodeMirror 6, uPlot, built with Vite) for the [GenSAS server](../doc/Server.md). See [doc/Web.md](../doc/Web.md).

```bash
cd web
npm ci
npm run dev        # http://localhost:5173, talks to a server at http://127.0.0.1:8765
npm test           # unit tests
npm run build      # static files in dist/
```

`npm run types` regenerates `src/api-types.ts` from `../api/openapi.yaml`; CI fails if the committed file is out of date. The end-to-end tests (`web/e2e`, Playwright and a real server) run with `python -m pytest web/e2e` after `npm run build`.
