import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";
import { useEffect, useRef } from "preact/hooks";
import { PALETTE, colorOf } from "./colors";
import type { ResultJson } from "./types";

export { PALETTE, colorOf };

interface Props {
  result: ResultJson;
  /** names to draw, in order; their color is the one of their position in `all` */
  visible: string[];
  all: string[];
  stacked: boolean;
}

function series(result: ResultJson, name: string): (number | null)[] {
  return result.variables[name] ?? result.discrete[name] ?? [];
}

function cssVar(name: string, fallback: string): string {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
}

/** Draws the times of the events as thin vertical lines. */
function eventLines(times: number[]): uPlot.Plugin {
  return {
    hooks: {
      draw: (u: uPlot) => {
        const ctx = u.ctx;
        ctx.save();
        ctx.strokeStyle = cssVar("--event-line", "rgba(120,120,120,0.6)");
        ctx.setLineDash([4, 4]);
        for (const t of times) {
          const x = Math.round(u.valToPos(t, "x", true));
          if (x < u.bbox.left || x > u.bbox.left + u.bbox.width) continue;
          ctx.beginPath();
          ctx.moveTo(x, u.bbox.top);
          ctx.lineTo(x, u.bbox.top + u.bbox.height);
          ctx.stroke();
        }
        ctx.restore();
      },
    },
  };
}

function build(host: HTMLElement, result: ResultJson, names: string[], all: string[], height: number, sync: string | null): uPlot {
  const axis = cssVar("--muted", "#666");
  const grid = cssVar("--grid", "#e5e5e5");
  const opts: uPlot.Options = {
    width: Math.max(host.clientWidth, 200),
    height,
    scales: { x: { time: false } },
    cursor: { sync: sync ? { key: sync } : undefined, drag: { x: true, y: false } },
    legend: { show: false },
    plugins: [eventLines(result.events.map((e) => e.time))],
    axes: [
      { stroke: axis, grid: { stroke: grid }, ticks: { stroke: grid }, label: "time" },
      { stroke: axis, grid: { stroke: grid }, ticks: { stroke: grid }, size: 60 },
    ],
    series: [
      {},
      ...names.map((n) => ({
        label: n,
        stroke: colorOf(all.indexOf(n)),
        width: 1.6,
        spanGaps: false,
        paths: result.discrete[n] ? uPlot.paths.stepped!({ align: 1 }) : undefined,
      })),
    ],
  };
  const data: uPlot.AlignedData = [result.time, ...names.map((n) => series(result, n))] as uPlot.AlignedData;
  return new uPlot(opts, data, host);
}

export function Plot({ result, visible, all, stacked }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const el = host.current!;
    el.replaceChildren();
    if (visible.length === 0) return;
    const plots: uPlot[] = [];
    if (stacked) {
      for (const name of visible) {
        const cell = document.createElement("div");
        cell.className = "plot-cell";
        const title = document.createElement("div");
        title.className = "plot-title";
        title.textContent = name;
        title.style.color = colorOf(all.indexOf(name));
        cell.append(title);
        el.append(cell);
        plots.push(build(cell, result, [name], all, 150, "gensas"));
      }
    } else {
      plots.push(build(el, result, visible, all, Math.max(260, Math.min(el.clientHeight || 380, 560)), null));
    }
    const resize = new ResizeObserver(() => {
      for (const p of plots) p.setSize({ width: Math.max(el.clientWidth, 200), height: p.height });
    });
    resize.observe(el);
    return () => {
      resize.disconnect();
      for (const p of plots) p.destroy();
    };
  }, [result, visible.join("\u0000"), all.join("\u0000"), stacked]);

  return <div class={stacked ? "plot stacked" : "plot"} ref={host} />;
}
