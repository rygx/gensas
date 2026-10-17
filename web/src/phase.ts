/** Geometry of the phase (x-y) and 3D plots; kept free of the DOM so that it can be tested. */

export type Series = (number | null)[];
export interface Range {
  min: number;
  max: number;
}

/** Smallest and largest finite value; a single value or no value gets a range around it. */
export function extent(series: Series[]): Range {
  let min = Infinity;
  let max = -Infinity;
  for (const s of series) {
    for (const v of s) {
      if (v === null || !Number.isFinite(v)) continue;
      if (v < min) min = v;
      if (v > max) max = v;
    }
  }
  if (min === Infinity) return { min: 0, max: 1 };
  if (min === max) {
    const pad = min === 0 ? 1 : Math.abs(min) * 0.1;
    return { min: min - pad, max: max + pad };
  }
  return { min, max };
}

/** Round tick positions inside [lo, hi], about `target` of them (steps of 1, 2 or 5 times a power of ten). */
export function niceTicks(lo: number, hi: number, target = 6): number[] {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / Math.max(1, target);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = (norm < 1.5 ? 1 : norm < 3.5 ? 2 : norm < 7.5 ? 5 : 10) * mag;
  const out: number[] = [];
  const first = Math.ceil(lo / step - 1e-9);
  for (let i = first; i * step <= hi + step * 1e-9; i++) {
    out.push(Number((i * step).toPrecision(12)));
  }
  return out;
}

export function formatTick(v: number): string {
  if (v === 0) return "0";
  const a = Math.abs(v);
  return a >= 1e5 || a < 1e-3 ? v.toExponential(1).replace("e+", "e") : String(Number(v.toPrecision(4)));
}

/** Maps v from [min, max] to [-1, 1]. */
export function normalize(v: number, r: Range): number {
  return ((v - r.min) / (r.max - r.min)) * 2 - 1;
}

export interface View3D {
  yaw: number;
  pitch: number;
  zoom: number;
}

export const DEFAULT_VIEW: View3D = { yaw: -0.65, pitch: 0.45, zoom: 1 };

export function clampView(v: View3D): View3D {
  return { yaw: v.yaw, pitch: Math.max(-Math.PI / 2, Math.min(Math.PI / 2, v.pitch)), zoom: Math.max(0.3, Math.min(6, v.zoom)) };
}

/**
 * Orthographic projection of a point of the cube [-1, 1]^3: turn about the vertical axis (z) by `yaw`, tilt
 * by `pitch`. Returns screen x to the right and y up, in the same units, and the depth (larger is nearer).
 */
export function project(p: [number, number, number], view: View3D): [number, number, number] {
  const cy = Math.cos(view.yaw);
  const sy = Math.sin(view.yaw);
  const cp = Math.cos(view.pitch);
  const sp = Math.sin(view.pitch);
  const x1 = p[0] * cy - p[1] * sy;
  const y1 = p[0] * sy + p[1] * cy;
  return [x1 * view.zoom, (p[2] * cp - y1 * sp) * view.zoom, (y1 * cp + p[2] * sp) * view.zoom];
}

/** The 12 edges of the cube [-1, 1]^3, as pairs of corners. */
export function cubeEdges(): [[number, number, number], [number, number, number]][] {
  const corners: [number, number, number][] = [];
  for (const x of [-1, 1]) for (const y of [-1, 1]) for (const z of [-1, 1]) corners.push([x, y, z]);
  const edges: [[number, number, number], [number, number, number]][] = [];
  for (let i = 0; i < corners.length; i++) {
    for (let j = i + 1; j < corners.length; j++) {
      const diff = corners[i].filter((v, k) => v !== corners[j][k]).length;
      if (diff === 1) edges.push([corners[i], corners[j]]);
    }
  }
  return edges;
}

/** Position 0..1 along a run of n samples as a color: from blue at the start to orange at the end. */
export function trailColor(t: number): string {
  const hue = 215 - 185 * Math.max(0, Math.min(1, t));
  return `hsl(${hue.toFixed(0)}, 75%, 48%)`;
}
