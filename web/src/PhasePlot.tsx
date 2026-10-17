import { useEffect, useRef, useState } from "preact/hooks";
import { colorOf } from "./colors";
import { DEFAULT_VIEW, type Range, type Series, type View3D, clampView, cubeEdges, extent, formatTick, niceTicks, normalize, project, trailColor } from "./phase";

interface Common {
  /** the values of every name, by name ("time" included) */
  data: Record<string, Series>;
  all: string[];
}

interface Props2D extends Common {
  x: string;
  ys: string[];
}

interface Props3D extends Common {
  x: string;
  y: string;
  z: string;
}

function cssVar(name: string, fallback: string): string {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
}

/** A canvas that follows the size of its container and the device pixel ratio; `draw` gets css pixels. */
function useCanvas(draw: (ctx: CanvasRenderingContext2D, w: number, h: number) => void, deps: unknown[]) {
  const canvas = useRef<HTMLCanvasElement>(null);
  useEffect(() => {
    const el = canvas.current!;
    const paint = () => {
      const w = Math.max(el.clientWidth, 200);
      const h = Math.max(el.clientHeight, 200);
      const ratio = window.devicePixelRatio || 1;
      el.width = Math.round(w * ratio);
      el.height = Math.round(h * ratio);
      const ctx = el.getContext("2d");
      if (!ctx) return;
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      ctx.clearRect(0, 0, w, h);
      draw(ctx, w, h);
    };
    paint();
    const observer = new ResizeObserver(paint);
    observer.observe(el);
    return () => observer.disconnect();
  }, deps);
  return canvas;
}

function polyline(ctx: CanvasRenderingContext2D, xs: Series, ys: Series, px: (v: number) => number, py: (v: number) => number) {
  ctx.beginPath();
  let pen = false;
  const n = Math.min(xs.length, ys.length);
  for (let i = 0; i < n; i++) {
    const x = xs[i];
    const y = ys[i];
    if (x === null || y === null || !Number.isFinite(x) || !Number.isFinite(y)) {
      pen = false;
      continue;
    }
    if (pen) ctx.lineTo(px(x), py(y));
    else ctx.moveTo(px(x), py(y));
    pen = true;
  }
  ctx.stroke();
}

/** Curves of several variables against one variable that is not necessarily the time. */
export function XYPlot({ data, all, x, ys }: Props2D) {
  const canvas = useCanvas(
    (ctx, w, h) => {
      const text = cssVar("--muted", "#666");
      const grid = cssVar("--grid", "#e5e5e5");
      const left = 64;
      const right = 14;
      const top = 12;
      const bottom = 40;
      const xr = extent([data[x] ?? []]);
      const yr = extent(ys.map((n) => data[n] ?? []));
      const px = (v: number) => left + ((v - xr.min) / (xr.max - xr.min)) * (w - left - right);
      const py = (v: number) => top + (1 - (v - yr.min) / (yr.max - yr.min)) * (h - top - bottom);
      ctx.font = "11px system-ui, sans-serif";
      ctx.fillStyle = text;
      ctx.strokeStyle = grid;
      ctx.lineWidth = 1;
      ctx.textAlign = "center";
      for (const t of niceTicks(xr.min, xr.max, Math.max(3, Math.floor((w - left) / 90)))) {
        ctx.beginPath();
        ctx.moveTo(px(t), top);
        ctx.lineTo(px(t), h - bottom);
        ctx.stroke();
        ctx.fillText(formatTick(t), px(t), h - bottom + 14);
      }
      ctx.textAlign = "right";
      for (const t of niceTicks(yr.min, yr.max, Math.max(3, Math.floor((h - top - bottom) / 40)))) {
        ctx.beginPath();
        ctx.moveTo(left, py(t));
        ctx.lineTo(w - right, py(t));
        ctx.stroke();
        ctx.fillText(formatTick(t), left - 6, py(t) + 4);
      }
      ctx.textAlign = "center";
      ctx.fillText(x, left + (w - left - right) / 2, h - 8);
      ctx.save();
      ctx.beginPath();
      ctx.rect(left, top, w - left - right, h - top - bottom);
      ctx.clip();
      ctx.lineWidth = 1.6;
      ctx.lineJoin = "round";
      for (const name of ys) {
        ctx.strokeStyle = colorOf(all.indexOf(name));
        polyline(ctx, data[x] ?? [], data[name] ?? [], px, py);
        // the start of the curve
        const x0 = (data[x] ?? []).find((v) => v !== null);
        const y0 = (data[name] ?? []).find((v) => v !== null);
        if (typeof x0 === "number" && typeof y0 === "number") {
          ctx.fillStyle = colorOf(all.indexOf(name));
          ctx.beginPath();
          ctx.arc(px(x0), py(y0), 3.5, 0, 2 * Math.PI);
          ctx.fill();
        }
      }
      ctx.restore();
    },
    [data, all.join("\u0000"), x, ys.join("\u0000")],
  );
  return <canvas ref={canvas} class="xy-canvas" role="img" aria-label={`${ys.join(", ")} against ${x}`} />;
}

/** One trajectory in three variables; drag to turn it, wheel to zoom, double click to reset. */
export function Plot3D({ data, x, y, z }: Props3D) {
  const [view, setView] = useState<View3D>(DEFAULT_VIEW);
  const drag = useRef<{ x: number; y: number } | null>(null);
  const names = [x, y, z];
  const ranges: Range[] = names.map((n) => extent([data[n] ?? []]));

  const canvas = useCanvas(
    (ctx, w, h) => {
      const text = cssVar("--muted", "#666");
      const grid = cssVar("--grid", "#cfd4da");
      const scale = Math.min(w, h) * 0.3;
      const cx = w / 2;
      const cy = h / 2;
      const at = (p: [number, number, number]): [number, number] => {
        const q = project(p, view);
        return [cx + q[0] * scale, cy - q[1] * scale];
      };
      ctx.lineWidth = 1;
      ctx.strokeStyle = grid;
      for (const [a, b] of cubeEdges()) {
        const [x1, y1] = at(a);
        const [x2, y2] = at(b);
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
      }
      // the three axes of the box, labeled at their far ends with the name and the range
      ctx.font = "11px system-ui, sans-serif";
      ctx.fillStyle = text;
      ctx.textAlign = "center";
      const axisEnds: [[number, number, number], [number, number, number]][] = [
        [[-1, -1, -1], [1, -1, -1]],
        [[-1, -1, -1], [-1, 1, -1]],
        [[-1, -1, -1], [-1, -1, 1]],
      ];
      axisEnds.forEach(([a, b], k) => {
        const [x1, y1] = at(a);
        const [x2, y2] = at(b);
        ctx.strokeStyle = text;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
        ctx.fillText(`${names[k]} ${formatTick(ranges[k].min)} … ${formatTick(ranges[k].max)}`, (x1 + x2) / 2 + (k === 2 ? -48 : 0), (y1 + y2) / 2 + (k === 2 ? 0 : 16));
      });
      const xs = data[x] ?? [];
      const ys = data[y] ?? [];
      const zs = data[z] ?? [];
      const n = Math.min(xs.length, ys.length, zs.length);
      const pts: ([number, number] | null)[] = [];
      for (let i = 0; i < n; i++) {
        const a = xs[i];
        const b = ys[i];
        const c = zs[i];
        if (a === null || b === null || c === null || !Number.isFinite(a) || !Number.isFinite(b) || !Number.isFinite(c)) pts.push(null);
        else pts.push(at([normalize(a, ranges[0]), normalize(b, ranges[1]), normalize(c, ranges[2])]));
      }
      // the color runs from the start to the end of the run, drawn in short pieces
      ctx.lineWidth = 1.4;
      ctx.lineJoin = "round";
      const pieces = Math.min(120, Math.max(1, n - 1));
      const per = Math.max(1, Math.ceil((n - 1) / pieces));
      for (let start = 0; start < n - 1; start += per) {
        ctx.strokeStyle = trailColor(start / Math.max(1, n - 1));
        ctx.beginPath();
        let pen = false;
        for (let i = start; i <= Math.min(n - 1, start + per); i++) {
          const p = pts[i];
          if (p === null) {
            pen = false;
            continue;
          }
          if (pen) ctx.lineTo(p[0], p[1]);
          else ctx.moveTo(p[0], p[1]);
          pen = true;
        }
        ctx.stroke();
      }
      const first = pts.find((p) => p !== null);
      if (first) {
        ctx.fillStyle = trailColor(0);
        ctx.beginPath();
        ctx.arc(first[0], first[1], 4, 0, 2 * Math.PI);
        ctx.fill();
      }
    },
    [data, x, y, z, view],
  );

  return (
    <div class="plot3d">
      <canvas
        ref={canvas}
        class="xy-canvas"
        role="img"
        aria-label={`${x}, ${y}, ${z} trajectory`}
        onPointerDown={(e) => {
          drag.current = { x: e.clientX, y: e.clientY };
          (e.target as HTMLElement).setPointerCapture?.(e.pointerId);
        }}
        onPointerMove={(e) => {
          if (!drag.current) return;
          const dx = e.clientX - drag.current.x;
          const dy = e.clientY - drag.current.y;
          drag.current = { x: e.clientX, y: e.clientY };
          setView((v) => clampView({ ...v, yaw: v.yaw + dx * 0.01, pitch: v.pitch + dy * 0.01 }));
        }}
        onPointerUp={() => (drag.current = null)}
        onPointerCancel={() => (drag.current = null)}
        onWheel={(e) => {
          e.preventDefault();
          setView((v) => clampView({ ...v, zoom: v.zoom * (e.deltaY < 0 ? 1.1 : 1 / 1.1) }));
        }}
        onDblClick={() => setView(DEFAULT_VIEW)}
      />
      <div class="muted plot3d-hint">drag to turn · wheel to zoom · double click to reset · the line goes from blue (start) to orange (end)</div>
    </div>
  );
}
