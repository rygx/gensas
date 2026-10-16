import type { Diagnostic, ModelInfo, Run } from "./types";

/** The text typed into a value field, kept as text so that "1e" or "-" can be typed on the way. */
export type Edits = Record<string, string>;

export function parseNumber(text: string): number | null {
  if (text.trim() === "") return null;
  const v = Number(text);
  return Number.isFinite(v) ? v : null;
}

/** The values to send: edited fields that parse as numbers. */
export function overrides(edits: Edits): Record<string, number> {
  const out: Record<string, number> = {};
  for (const [k, v] of Object.entries(edits)) {
    const n = parseNumber(v);
    if (n !== null) out[k] = n;
  }
  return out;
}

interface ParamsProps {
  info: ModelInfo | null;
  params: Edits;
  start: Edits;
  setParams: (e: Edits) => void;
  setStart: (e: Edits) => void;
}

function EditRow(props: { name: string; shown: string; hint: string; edit: string | undefined; disabled?: boolean; onEdit: (v: string | null) => void; unit: string; text: string }) {
  const changed = props.edit !== undefined;
  const invalid = changed && parseNumber(props.edit!) === null;
  return (
    <tr class={changed ? "changed" : ""}>
      <td class="name" title={props.text}>
        {props.name}
      </td>
      <td>
        <input
          type="text"
          inputMode="decimal"
          value={changed ? props.edit : props.shown}
          disabled={props.disabled}
          class={invalid ? "invalid" : ""}
          aria-label={props.name}
          title={props.hint}
          onInput={(e) => props.onEdit((e.target as HTMLInputElement).value === props.shown ? null : (e.target as HTMLInputElement).value)}
        />
      </td>
      <td class="unit">{props.unit}</td>
      <td class="reset">
        {changed && (
          <button class="link" title="back to the value in the model" onClick={() => props.onEdit(null)}>
            reset
          </button>
        )}
      </td>
    </tr>
  );
}

export function ParamsPanel({ info, params, start, setParams, setStart }: ParamsProps) {
  if (!info) return <p class="muted">Parameters and start values appear when the model has no errors.</p>;
  const shownParams = info.parameters;
  const states = info.variables.filter((v) => v.kind === "state" && !v.input);
  const set = (edits: Edits, update: (e: Edits) => void) => (name: string) => (v: string | null) => {
    const next = { ...edits };
    if (v === null) delete next[name];
    else next[name] = v;
    update(next);
  };
  const anyEdit = Object.keys(params).length + Object.keys(start).length > 0;
  return (
    <div class="params">
      {shownParams.length === 0 && states.length === 0 && <p class="muted">The model has no parameters or states.</p>}
      {shownParams.length > 0 && (
        <table>
          <caption>Parameters</caption>
          <tbody>
            {shownParams.map((p) => (
              <EditRow
                key={p.name}
                name={p.name}
                shown={String(p.value)}
                hint={p.constant ? "constants cannot be changed" : "value in the model: " + p.value}
                edit={params[p.name]}
                disabled={p.constant}
                onEdit={set(params, setParams)(p.name)}
                unit={p.unit}
                text={p.description}
              />
            ))}
          </tbody>
        </table>
      )}
      {states.length > 0 && (
        <table>
          <caption>Start values of the states</caption>
          <tbody>
            {states.map((v) => (
              <EditRow
                key={v.name}
                name={v.name}
                shown={v.start === null ? "" : String(v.start)}
                hint={v.start === null ? "no start value in the model" : "start value in the model: " + v.start}
                edit={start[v.name]}
                onEdit={set(start, setStart)(v.name)}
                unit={v.unit}
                text={v.description}
              />
            ))}
          </tbody>
        </table>
      )}
      {anyEdit && (
        <button
          class="link"
          onClick={() => {
            setParams({});
            setStart({});
          }}
        >
          reset all values
        </button>
      )}
    </div>
  );
}

export function DiagnosticList({ items, onJump }: { items: Diagnostic[]; onJump: (line: number) => void }) {
  if (items.length === 0) return null;
  return (
    <ul class="diagnostics">
      {items.map((d, i) => (
        <li key={i} class={d.severity}>
          <button class="link" disabled={d.line < 1} onClick={() => onJump(d.line)}>
            {d.line > 0 ? `line ${d.line}` : "model"}
          </button>
          <span>{d.message}</span>
        </li>
      ))}
    </ul>
  );
}

const STATE_LABEL: Record<Run["state"], string> = {
  queued: "waiting for a free worker",
  running: "running",
  completed: "completed",
  terminated: "stopped by terminate()",
  failed: "failed",
  cancelled: "cancelled",
};

export function RunStatus({ run }: { run: Run | null }) {
  if (!run) return null;
  const fraction = run.progress.end > 0 ? Math.min(1, run.progress.time / run.progress.end) : 0;
  const active = run.state === "queued" || run.state === "running";
  return (
    <div class={"run-status " + run.state}>
      <div class="line">
        <strong>{STATE_LABEL[run.state]}</strong>
        {run.elapsedSeconds !== null && run.elapsedSeconds !== undefined && <span class="muted"> {run.elapsedSeconds.toFixed(2)} s</span>}
        {run.termination && <span class="muted"> at t = {run.termination.time.toPrecision(5)}: {run.termination.message}</span>}
      </div>
      {active && <progress max={1} value={fraction} aria-label="progress" />}
      {run.error && <div class="error-text">{run.error}</div>}
      {run.warnings.length > 0 && (
        <ul class="warnings">
          {run.warnings.map((w, i) => (
            <li key={i}>{w}</li>
          ))}
        </ul>
      )}
    </div>
  );
}
