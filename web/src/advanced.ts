function parseNumber(text: string): number | null {
  if (text.trim() === "") return null;
  const v = Number(text);
  return Number.isFinite(v) ? v : null;
}

/**
 * Solver settings of the advanced panel. The ranges are the ones the engine accepts (doc/JobFormat.md);
 * the engine moves a value outside them to the nearest limit with a warning, the client does not send it.
 */
export interface Field {
  key: "segment" | "level" | "alphaTol" | "errorTol" | "outputStep" | "smoothEps";
  label: string;
  help: string;
  min: number;
  max: number;
  /** the engine's default, shown as the placeholder; null when it depends on the other settings */
  fallback: string;
  integer?: boolean;
}

export const FIELDS: Field[] = [
  { key: "segment", label: "Longest step", help: "Largest length of one solver segment (s). Shorter steps follow fast changes more closely and cost more time.", min: 1e-4, max: 10, fallback: "1" },
  { key: "level", label: "Series order", help: "Order of the series solution in each step. A higher order allows longer steps for smooth models.", min: 1, max: 1000, fallback: "30", integer: true },
  { key: "alphaTol", label: "Step shrink tolerance", help: "How finely the step end is searched when the error is too large. At most the longest step.", min: 1e-9, max: 10, fallback: "0.001" },
  { key: "errorTol", label: "Error tolerance", help: "Largest residual of the equations accepted in a step. Smaller is more accurate and slower.", min: 1e-10, max: 1e-2, fallback: "1e-5" },
  { key: "outputStep", label: "Output spacing", help: "Spacing of the result samples (s). The solver steps are not affected. Default: one thousandth of the end time.", min: 1e-3, max: 1, fallback: "automatic" },
  { key: "smoothEps", label: "Smoothing width", help: "Width of the rounded corner when abs, min, max and sign are smoothed.", min: Number.MIN_VALUE, max: Infinity, fallback: "1e-6" },
];

/** What was typed, by key; "piecewise" holds "events" or "smooth". An absent key means the default. */
export type AdvancedEdits = Record<string, string>;

export type Piecewise = "events" | "smooth";

export function piecewiseOf(edits: AdvancedEdits): Piecewise {
  return edits.piecewise === "smooth" ? "smooth" : "events";
}

/** Problems by key for the typed values. */
export function advancedErrors(edits: AdvancedEdits): Record<string, string> {
  const out: Record<string, string> = {};
  const segment = parseNumber(edits.segment ?? "") ?? 1;
  for (const f of FIELDS) {
    const text = edits[f.key];
    if (text === undefined || text.trim() === "") continue;
    if (f.key === "smoothEps" && piecewiseOf(edits) !== "smooth") continue;
    const v = parseNumber(text);
    const max = f.key === "alphaTol" ? Math.min(f.max, segment) : f.max;
    if (v === null) out[f.key] = "not a number";
    else if (f.integer && !Number.isInteger(v)) out[f.key] = "must be a whole number";
    else if (f.key === "smoothEps" ? v <= 0 : v < f.min || v > max) out[f.key] = f.key === "smoothEps" ? "must be positive" : `between ${f.min} and ${max}`;
  }
  return out;
}

export interface AdvancedRequest {
  solver: Record<string, number>;
  model: { piecewise?: Piecewise; smoothEps?: number };
  outputStep?: number;
}

/** The fields of the run request that the typed values change. */
export function advancedRequest(edits: AdvancedEdits): AdvancedRequest {
  const errors = advancedErrors(edits);
  const num = (key: string) => (errors[key] === undefined ? parseNumber(edits[key] ?? "") : null);
  const solver: Record<string, number> = {};
  for (const key of ["segment", "level", "alphaTol", "errorTol"]) {
    const v = num(key);
    if (v !== null) solver[key] = v;
  }
  const model: AdvancedRequest["model"] = {};
  if (piecewiseOf(edits) === "smooth") {
    model.piecewise = "smooth";
    const eps = num("smoothEps");
    if (eps !== null) model.smoothEps = eps;
  }
  const outputStep = num("outputStep");
  return { solver, model, ...(outputStep !== null ? { outputStep } : {}) };
}

export function isAdvancedChanged(edits: AdvancedEdits): boolean {
  return Object.keys(edits).length > 0;
}
