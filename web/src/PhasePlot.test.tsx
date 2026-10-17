import { render } from "preact";
import { act } from "preact/test-utils";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { Plot3D, XYPlot } from "./PhasePlot";

// A canvas context that only counts what is drawn.
function fakeContext() {
  const calls: Record<string, number> = {};
  const texts: string[] = [];
  const ctx = new Proxy(
    {},
    {
      get(_t, prop: string) {
        if (prop === "fillText") return (t: string) => texts.push(t);
        return (..._a: unknown[]) => {
          calls[prop] = (calls[prop] ?? 0) + 1;
        };
      },
      set() {
        return true;
      },
    },
  );
  return { ctx, calls, texts };
}

let host: HTMLElement;
let fake: ReturnType<typeof fakeContext>;
beforeEach(() => {
  fake = fakeContext();
  vi.spyOn(HTMLCanvasElement.prototype, "getContext").mockImplementation((() => fake.ctx) as never);
  vi.stubGlobal("ResizeObserver", class { observe() {} disconnect() {} });
  host = document.createElement("div");
  document.body.append(host);
});
afterEach(() => {
  act(() => render(null, host));
  host.remove();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

const data = { time: [0, 1, 2, 3], x: [0, 1, null, 0], y: [1, 0, -1, 0], z: [0, 0.5, 1, 0.5] };

describe("XYPlot", () => {
  it("draws a curve per variable and names the x axis", () => {
    act(() => render(<XYPlot data={data} all={["x", "y", "z"]} x="x" ys={["y", "z"]} />, host));
    expect(host.querySelector("canvas")!.getAttribute("aria-label")).toBe("y, z against x");
    expect(fake.calls.stroke).toBeGreaterThan(4); // grid lines and two curves
    expect(fake.calls.arc).toBe(2); // the starts
    expect(fake.texts).toContain("x");
  });

  it("draws nothing for unknown names without failing", () => {
    act(() => render(<XYPlot data={data} all={["x"]} x="nope" ys={["also-nope"]} />, host));
    expect(host.querySelector("canvas")).not.toBeNull();
  });
});

describe("Plot3D", () => {
  it("draws the box, the three axes and the trajectory", () => {
    act(() => render(<Plot3D data={data} all={["x", "y", "z"]} x="x" y="y" z="z" />, host));
    expect(fake.calls.stroke).toBeGreaterThan(15);
    expect(fake.texts.some((t) => t.startsWith("x "))).toBe(true);
    expect(fake.texts.some((t) => t.startsWith("z "))).toBe(true);
    expect(fake.calls.arc).toBe(1);
  });

  it("turns with the pointer, zooms with the wheel and resets on a double click", () => {
    act(() => render(<Plot3D data={data} all={["x", "y", "z"]} x="x" y="y" z="z" />, host));
    const canvas = host.querySelector("canvas")!;
    const strokes = () => fake.calls.stroke;
    const before = strokes();
    act(() => {
      canvas.dispatchEvent(new PointerEvent("pointerdown", { clientX: 10, clientY: 10, bubbles: true }));
      canvas.dispatchEvent(new PointerEvent("pointermove", { clientX: 50, clientY: 30, bubbles: true }));
      canvas.dispatchEvent(new PointerEvent("pointerup", { bubbles: true }));
    });
    expect(strokes()).toBeGreaterThan(before); // redrawn
    const afterTurn = strokes();
    act(() => {
      canvas.dispatchEvent(new WheelEvent("wheel", { deltaY: -100, bubbles: true, cancelable: true }));
    });
    expect(strokes()).toBeGreaterThan(afterTurn);
    const afterZoom = strokes();
    act(() => {
      canvas.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
    });
    expect(strokes()).toBeGreaterThan(afterZoom);
  });
});
