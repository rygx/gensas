# Web client

`web/` is a browser client for the [GenSAS server](Server.md): edit a Modelica model, see errors with line numbers while typing, change parameters and start values, run, cancel, and plot the curves. It is a static site (no server code of its own), so it can be hosted anywhere, and it talks to a server on your machine or on a remote one through the [API](../api/openapi.yaml).

## Using it

1. Start a server (see [Server](Server.md)): `python -m gensas_server --static-dir ../web/dist` serves the built client at the server's own address, so one command gives you the whole application at `http://127.0.0.1:8764/`.
2. Or serve `web/dist` from any web server and open **Server** in the header to enter the address (and the access token, if the server wants one). The choice is remembered in the browser.
3. Pick an example or type a model. After a short pause the model is checked: errors and warnings appear in the margin and in the list under the editor (click `line N` to jump there). **Run** (or Ctrl+Enter) starts a simulation; **Cancel** stops it.
4. **Parameters and start values** list what the model declares. Type a new value to change it for the next run (the field is highlighted, **reset** puts the model's value back); the model text is not touched.

   **Advanced solver settings** (below the parameters) change how the solver works for the next run: the longest step, the order of the series, the step shrink tolerance, the error tolerance, the spacing of the result samples, and whether `abs`, `min`, `max` and `sign` are exact (events) or smoothed (with a smoothing width). An empty field uses the default. Each field has the range the engine accepts; a value outside it is marked and the Run button stays disabled. Smaller tolerances and shorter steps are more accurate and slower.
5. The plot shows the variables ticked above it. **One plot per variable** stacks them when their scales differ. Dashed vertical lines are events. Drag in the plot to zoom, double-click to reset. **CSV** and **MAT** download the full result.

## Development

```bash
cd web
npm ci
npm run dev       # hot-reloading page on :5173; allowed by the server's default CORS rule for local origins
npm test          # unit tests (API client, event stream parser)
npm run build
```

The client is outside Bazel. `src/api-types.ts` is generated from `api/openapi.yaml` (`npm run types`); CI regenerates it and fails on a difference. The Web client workflow also runs end-to-end tests in a real Chromium against a real server (`web/e2e`).

Limits: one model in one file, one result at a time, and the plot shows at most about 4000 samples per variable (the server thins larger results and keeps the extremes and event times).
