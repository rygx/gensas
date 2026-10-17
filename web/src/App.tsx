import { useCallback, useEffect, useMemo, useRef, useState } from "preact/hooks";
import { Api, ApiError } from "./api";
import { Editor } from "./Editor";
import { type AdvancedEdits, advancedErrors, advancedRequest, isAdvancedChanged } from "./advanced";
import { AdvancedPanel, DiagnosticList, type Edits, ParamsPanel, RunStatus, overrides } from "./Panels";
import { Plot, colorOf } from "./Plot";
import { Plot3D, XYPlot } from "./PhasePlot";
import { type Settings, loadSettings, loadSource, saveSettings, saveSource } from "./settings";
import { type CheckResponse, type Diagnostic, type Example, FINISHED, HAS_RESULT, type ModelInfo, type ResultJson, type Run, type Version } from "./types";

const SUPPORTED_API = 1;
const MAX_PLOT_POINTS = 4000;

function describe(err: unknown): string {
  if (err instanceof ApiError) return err.message;
  return err instanceof Error ? err.message : String(err);
}

function download(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function App() {
  const saved = useMemo(loadSource, []);
  const [settings, setSettings] = useState<Settings>(loadSettings);
  const [showSettings, setShowSettings] = useState(false);
  const [version, setVersion] = useState<Version | null>(null);
  const [connError, setConnError] = useState<string | null>(null);
  const [examples, setExamples] = useState<Example[]>([]);

  const [source, setSource] = useState(saved?.source ?? "");
  const [endTime, setEndTime] = useState(String(saved?.endTime ?? 5));
  const [checked, setChecked] = useState<CheckResponse | null>(null);
  const [info, setInfo] = useState<ModelInfo | null>(null);
  const [models, setModels] = useState<string[]>([]);
  const [className, setClassName] = useState("");
  const [paramEdits, setParamEdits] = useState<Edits>({});
  const [startEdits, setStartEdits] = useState<Edits>({});
  const [advanced, setAdvanced] = useState<AdvancedEdits>({});

  const [run, setRun] = useState<Run | null>(null);
  const [result, setResult] = useState<ResultJson | null>(null);
  const [runError, setRunError] = useState<string | null>(null);
  const [selected, setSelected] = useState<string[]>([]);
  const [stacked, setStacked] = useState(false);
  const [view, setView] = useState<"time" | "xy" | "3d">("time");
  const [axes, setAxes] = useState({ x: "time", a: "", b: "", c: "" });
  const [jump, setJump] = useState<{ line: number } | null>(null);
  const abort = useRef<AbortController | null>(null);

  const api = useMemo(() => new Api(settings.serverUrl, settings.token), [settings]);

  // ---- connection ----
  useEffect(() => {
    let current = true;
    setVersion(null);
    setConnError(null);
    (async () => {
      try {
        const v = await api.version();
        const list = await api.examples();
        if (!current) return;
        setVersion(v);
        setExamples(list);
        if (!loadSource() && list.length > 0) {
          setSource(list[0].source);
          setEndTime(String(list[0].endTime ?? 5));
        }
      } catch (err) {
        if (current) setConnError(describe(err));
      }
    })();
    return () => {
      current = false;
    };
  }, [api]);

  // ---- checking while typing ----
  useEffect(() => {
    if (!version || source.trim() === "") return;
    let current = true;
    const timer = setTimeout(async () => {
      try {
        const r = await api.check({ source, className: className || undefined });
        if (!current) return;
        setChecked(r);
        if (r.info) {
          setModels(r.info.models ?? []);
          if (r.info.complete !== false) setInfo(r.info);
        }
      } catch (err) {
        if (current) setChecked({ ok: false, diagnostics: [{ severity: "error", line: 0, column: 0, message: describe(err), source: "server" }], info: null });
      }
    }, 500);
    return () => {
      current = false;
      clearTimeout(timer);
    };
  }, [source, className, version, api]);

  // a model that the text no longer has is dropped
  useEffect(() => {
    if (className !== "" && models.length > 0 && !models.includes(className)) setClassName("");
  }, [models, className]);

  useEffect(() => saveSource(source, Number(endTime) || 5), [source, endTime]);

  // edits of names the model no longer has are dropped
  useEffect(() => {
    if (!info) return;
    const keep = (edits: Edits, names: Set<string>) => {
      const next = Object.fromEntries(Object.entries(edits).filter(([k]) => names.has(k)));
      return Object.keys(next).length === Object.keys(edits).length ? edits : next;
    };
    setParamEdits((e) => keep(e, new Set(info.parameters.map((p) => p.name))));
    setStartEdits((e) => keep(e, new Set(info.variables.map((v) => v.name))));
  }, [info]);

  const diagnostics: Diagnostic[] = checked?.diagnostics ?? [];
  const errors = diagnostics.filter((d) => d.severity === "error");
  const endValue = Number(endTime);
  const endValid = Number.isFinite(endValue) && endValue > 0;
  const advancedValid = Object.keys(advancedErrors(advanced)).length === 0;
  const active = run !== null && (run.state === "queued" || run.state === "running");
  const canRun = version !== null && !active && endValid && advancedValid && source.trim() !== "" && (checked === null || checked.ok);

  const start = useCallback(async () => {
    setRunError(null);
    setResult(null);
    abort.current?.abort();
    const controller = new AbortController();
    abort.current = controller;
    try {
      const adv = advancedRequest(advanced);
      const created = await api.createRun({
        model: { source, className: className || undefined, parameters: overrides(paramEdits), start: overrides(startEdits), ...adv.model },
        solver: { ...adv.solver, endTime: endValue },
        output: { step: adv.outputStep ?? Math.min(1, Math.max(0.001, endValue / 1000)) },
      });
      setRun(created);
      const final = await api.follow(created.id, setRun, controller.signal);
      setRun(final);
      if (HAS_RESULT.includes(final.state)) {
        const r = await api.result(final.id, { maxPoints: MAX_PLOT_POINTS });
        setResult(r);
        const names = [...Object.keys(r.variables), ...Object.keys(r.discrete)];
        setSelected(names.length <= 6 ? names : names.slice(0, 4));
        // the variables of the phase portrait and the 3D plot: the first ones
        setAxes({ x: names.length > 1 ? names[0] : "time", a: names[0] ?? "", b: names[1] ?? names[0] ?? "", c: names[2] ?? names[1] ?? names[0] ?? "" });
      }
    } catch (err) {
      if ((err as Error).name !== "AbortError") setRunError(describe(err));
    }
  }, [api, source, className, paramEdits, startEdits, advanced, endValue]);

  const cancel = useCallback(async () => {
    if (run) {
      try {
        await api.cancelRun(run.id);
      } catch (err) {
        setRunError(describe(err));
      }
    }
  }, [api, run]);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.ctrlKey || e.metaKey) && e.key === "Enter") {
        e.preventDefault();
        if (canRun) start();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [canRun, start]);

  const allNames = result ? [...Object.keys(result.variables), ...Object.keys(result.discrete)] : [];
  const data = useMemo(() => (result ? ({ time: result.time, ...result.variables, ...result.discrete } as Record<string, (number | null)[]>) : {}), [result]);
  const toggle = (name: string) => setSelected((s) => (s.includes(name) ? s.filter((n) => n !== name) : [...s, name]));
  const saveConnection = (s: Settings) => {
    saveSettings(s);
    setSettings(s);
    setShowSettings(false);
  };

  return (
    <div class="app">
      <header>
        <h1>GenSAS</h1>
        <select
          aria-label="example"
          value=""
          onChange={(e) => {
            const ex = examples.find((x) => x.id === (e.target as HTMLSelectElement).value);
            if (ex) {
              setSource(ex.source);
              setEndTime(String(ex.endTime ?? 5));
              setClassName("");
              setParamEdits({});
              setStartEdits({});
              setResult(null);
              setRun(null);
              setChecked(null);
            }
          }}
          disabled={examples.length === 0}
        >
          <option value="">Load an example…</option>
          {examples.map((e) => (
            <option key={e.id} value={e.id} title={e.description}>
              {e.title}
            </option>
          ))}
        </select>
        <span class="spacer" />
        <span class={"conn " + (version ? "ok" : connError ? "bad" : "")} title={connError ?? ""}>
          {version ? `server ${version.server} · engine ${version.engine.version}` : connError ? "not connected" : "connecting…"}
        </span>
        <button onClick={() => setShowSettings(!showSettings)} aria-expanded={showSettings}>
          Server
        </button>
      </header>

      {showSettings && <ConnectionForm settings={settings} onSave={saveConnection} required={version?.authRequired ?? false} />}
      {connError && (
        <div class="banner bad">
          {connError}.{" "}
          {connError.includes("token") || connError.includes("bearer")
            ? "Enter the access token under Server."
            : "Start the server with `python -m gensas_server` or change the address under Server."}
        </div>
      )}
      {version && version.api !== SUPPORTED_API && (
        <div class="banner warn">This server speaks API version {version.api}; this page expects {SUPPORTED_API}. Some things may not work.</div>
      )}

      <main>
        <section class="left">
          <Editor value={source} onChange={setSource} diagnostics={diagnostics} jumpTo={jump} />
          <div class="status-line">
            {models.length > 1 && (
              <select aria-label="model" value={className || models[models.length - 1]} onChange={(e) => setClassName((e.target as HTMLSelectElement).value)}>
                {models.map((m) => (
                  <option key={m} value={m}>
                    {m}
                  </option>
                ))}
              </select>
            )}
            {checked === null ? (
              <span class="muted">{version ? "checking…" : ""}</span>
            ) : checked.ok ? (
              <span class="ok-text">
                ✓ {info?.name ?? "model"} {info?.hybrid ? "(with events)" : ""}
              </span>
            ) : (
              <span class="error-text">
                {errors.length} error{errors.length === 1 ? "" : "s"}
              </span>
            )}
          </div>
          <DiagnosticList items={diagnostics} onJump={(line) => setJump({ line })} />
        </section>

        <section class="right">
          <div class="run-bar">
            <label>
              End time
              <input type="text" inputMode="decimal" class={endValid ? "" : "invalid"} value={endTime} onInput={(e) => setEndTime((e.target as HTMLInputElement).value)} aria-label="end time" />
            </label>
            {active ? (
              <button class="danger" onClick={cancel}>
                Cancel
              </button>
            ) : (
              <button class="primary" onClick={start} disabled={!canRun} title="Ctrl+Enter">
                Run
              </button>
            )}
          </div>
          {runError && <div class="banner bad">{runError}</div>}
          <RunStatus run={run} />
          {run && run.diagnostics.length > 0 && FINISHED.includes(run.state) && run.state === "failed" && <DiagnosticList items={run.diagnostics} onJump={(line) => setJump({ line })} />}

          <details open>
            <summary>Parameters and start values</summary>
            <ParamsPanel info={info} params={paramEdits} start={startEdits} setParams={setParamEdits} setStart={setStartEdits} />
          </details>

          <details>
            <summary>Advanced solver settings{isAdvancedChanged(advanced) ? " (changed)" : ""}</summary>
            <AdvancedPanel edits={advanced} setEdits={setAdvanced} />
          </details>

          {result && run && (
            <div class="results">
              <div class="result-bar">
                <strong>Result</strong>
                <select aria-label="plot type" value={view} onChange={(e) => setView((e.target as HTMLSelectElement).value as typeof view)}>
                  <option value="time">against time</option>
                  <option value="xy">against a variable (x–y)</option>
                  <option value="3d">3D trajectory</option>
                </select>
                {view === "time" && (
                  <label>
                    <input type="checkbox" checked={stacked} onChange={(e) => setStacked((e.target as HTMLInputElement).checked)} /> one plot per variable
                  </label>
                )}
                {view === "xy" && (
                  <label>
                    x axis
                    <AxisSelect label="x axis" value={axes.x} names={["time", ...allNames]} onChange={(v) => setAxes({ ...axes, x: v })} />
                  </label>
                )}
                {view === "3d" && (
                  <>
                    <label>
                      x
                      <AxisSelect label="x variable" value={axes.a} names={allNames} onChange={(v) => setAxes({ ...axes, a: v })} />
                    </label>
                    <label>
                      y
                      <AxisSelect label="y variable" value={axes.b} names={allNames} onChange={(v) => setAxes({ ...axes, b: v })} />
                    </label>
                    <label>
                      z
                      <AxisSelect label="z variable" value={axes.c} names={allNames} onChange={(v) => setAxes({ ...axes, c: v })} />
                    </label>
                  </>
                )}
                <span class="spacer" />
                <button class="link" onClick={async () => download(await api.resultFile(run.id, "csv"), `${info?.name ?? "result"}.csv`)}>
                  CSV
                </button>
                <button class="link" onClick={async () => download(await api.resultFile(run.id, "mat"), `${info?.name ?? "result"}.mat`)}>
                  MAT
                </button>
              </div>
              {view !== "3d" && (
              <ul class="legend">
                {allNames.map((n, i) => (
                  <li key={n}>
                    <label>
                      <input type="checkbox" checked={selected.includes(n)} onChange={() => toggle(n)} />
                      <span class="swatch" style={{ background: colorOf(i) }} />
                      {n}
                    </label>
                  </li>
                ))}
              </ul>
              )}
              {view === "time" && <Plot result={result} visible={allNames.filter((n) => selected.includes(n))} all={allNames} stacked={stacked} />}
              {view === "xy" && <XYPlot data={data} all={allNames} x={axes.x} ys={allNames.filter((n) => selected.includes(n) && n !== axes.x)} />}
              {view === "3d" && <Plot3D data={data} all={allNames} x={axes.a} y={axes.b} z={axes.c} />}
              {result.events.length > 0 && (
                <details>
                  <summary>{result.events.length} event{result.events.length === 1 ? "" : "s"}</summary>
                  <table class="events">
                    <thead>
                      <tr>
                        <th>time</th>
                        <th>changes</th>
                      </tr>
                    </thead>
                    <tbody>
                      {result.events.slice(0, 200).map((e, i) => (
                        <tr key={i}>
                          <td>{e.time.toPrecision(6)}</td>
                          <td>
                            {[...Object.entries(e.assignments).map(([k, v]) => `${k} := ${v}`), ...Object.entries(e.reinits).map(([k, v]) => `reinit ${k} = ${v}`)].join(", ")}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </details>
              )}
            </div>
          )}
        </section>
      </main>
    </div>
  );
}

function AxisSelect(props: { label: string; value: string; names: string[]; onChange: (v: string) => void }) {
  return (
    <select aria-label={props.label} value={props.value} onChange={(e) => props.onChange((e.target as HTMLSelectElement).value)}>
      {props.names.map((n) => (
        <option key={n} value={n}>
          {n}
        </option>
      ))}
    </select>
  );
}

function ConnectionForm(props: { settings: Settings; onSave: (s: Settings) => void; required: boolean }) {
  const [url, setUrl] = useState(props.settings.serverUrl);
  const [token, setToken] = useState(props.settings.token);
  return (
    <form
      class="connection"
      onSubmit={(e) => {
        e.preventDefault();
        props.onSave({ serverUrl: url.trim(), token });
      }}
    >
      <label>
        Server address
        <input type="url" value={url} onInput={(e) => setUrl((e.target as HTMLInputElement).value)} required />
      </label>
      <label>
        Access token {props.required ? "(required)" : "(if the server asks for one)"}
        <input type="password" value={token} onInput={(e) => setToken((e.target as HTMLInputElement).value)} autocomplete="off" />
      </label>
      <button class="primary" type="submit">
        Connect
      </button>
    </form>
  );
}
