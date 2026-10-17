import { act } from "preact/test-utils";
import { render } from "preact";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

// CodeMirror and uPlot need a real layout engine and a canvas; the end-to-end tests cover them.
vi.mock("./Editor", () => ({
  Editor: (p: { value: string; onChange: (t: string) => void }) => <textarea data-testid="editor" value={p.value} onInput={(e) => p.onChange((e.target as HTMLTextAreaElement).value)} />,
}));
vi.mock("./Plot", () => ({
  colorOf: () => "#000",
  Plot: (p: { visible: string[]; stacked: boolean }) => <div data-testid="plot">{p.visible.join(",")}{p.stacked ? " stacked" : ""}</div>,
}));

vi.mock("./PhasePlot", () => ({
  XYPlot: (p: { x: string; ys: string[] }) => <div data-testid="xy">{p.ys.join(",")} against {p.x}</div>,
  Plot3D: (p: { x: string; y: string; z: string }) => <div data-testid="p3d">{[p.x, p.y, p.z].join(",")}</div>,
}));

import { App } from "./App";

const version = { server: "0.1.0", api: 1, engine: { version: "GenSAS", apiVersion: 2 }, authRequired: false, limits: { maxSourceBytes: 1, maxConcurrentRuns: 1, runTimeoutSeconds: 1 } };
const examples = [
  { id: "decay", title: "Decay", description: "", source: "model Decay\nend Decay;\n", endTime: 4 },
  { id: "two", title: "Two", description: "", source: "model Two\nend Two;\n", endTime: 2 },
];
const info = {
  name: "Decay",
  description: "",
  hybrid: false,
  parameters: [{ name: "k", type: "Real", unit: "", description: "", value: 1, constant: false, line: 2 }],
  variables: [{ name: "x", kind: "state", type: "Real", unit: "", description: "", start: 1, input: false, line: 3 }],
};
const run = (state: string, extra = {}) => ({ id: "r1", kind: "modelica", state, createdAt: 0, startedAt: 0, finishedAt: null, progress: { time: 0, end: 4 }, error: null, diagnostics: [], warnings: [], summary: null, ...extra });
const result = { time: [0, 1], variables: { x: [1, 0.5], y: [0, 1] }, discrete: {}, events: [{ time: 0.5, assignments: { d: 1 }, reinits: {} }], terminated: null };

let calls: { method: string; path: string; body?: any }[] = [];
let checkReply: any;
let runFinal: any;

function json(body: unknown, status = 200) {
  return new Response(JSON.stringify(body), { status });
}

beforeEach(() => {
  localStorage.clear();
  calls = [];
  checkReply = { ok: true, diagnostics: [], info };
  runFinal = run("completed", { summary: { samples: 2, variables: ["x", "y"], discrete: [], events: 1, endTime: 1 } });
  vi.spyOn(globalThis, "fetch").mockImplementation(async (input, init) => {
    const url = new URL(String(input));
    const method = init?.method ?? "GET";
    const body = init?.body ? JSON.parse(String(init.body)) : undefined;
    calls.push({ method, path: url.pathname + url.search, body });
    if (url.pathname === "/version") return json(version);
    if (url.pathname === "/examples") return json(examples);
    if (url.pathname === "/check") return json(checkReply);
    if (url.pathname === "/runs" && method === "POST") return json(run("queued"), 202);
    if (url.pathname === "/runs/r1/events") return new Response(`event: state\ndata: ${JSON.stringify(run("running"))}\n\nevent: done\ndata: ${JSON.stringify(runFinal)}\n\n`);
    if (url.pathname === "/runs/r1/result") return json(result);
    if (url.pathname === "/runs/r1" && method === "DELETE") return json(run("cancelled"));
    return json({ error: { code: "not_found", message: "no" } }, 404);
  });
});

let host: HTMLElement;
afterEach(async () => {
  await act(async () => render(null, host)); // runs the cleanups of the effects before the environment goes away
  host.remove();
  vi.restoreAllMocks();
});

async function mount(): Promise<HTMLElement> {
  host = document.createElement("div");
  document.body.append(host);
  await act(async () => {
    render(<App />, host);
  });
  await settle();
  return host;
}

// Waits in small steps: preact's act() only renders when it ends, so one long act would start the effects
// (and the timers of the debounced check) too late.
async function settle(ms = 1100) {
  for (let t = 0; t < ms; t += 50) {
    await act(async () => {
      await new Promise((r) => setTimeout(r, 50));
    });
  }
}

const click = (el: Element | null) => act(() => (el as HTMLElement).click());
const button = (name: string) => Array.from(host.querySelectorAll("button")).find((b) => b.textContent?.trim() === name) ?? null;

describe("App", () => {
  it("connects, shows the first example and its check", async () => {
    const el = await mount();
    expect(el.querySelector(".conn.ok")!.textContent).toContain("server 0.1.0");
    expect((el.querySelector("[data-testid=editor]") as HTMLTextAreaElement).value).toContain("model Decay");
    expect(el.querySelector(".status-line .ok-text")!.textContent).toContain("Decay");
    expect(el.querySelectorAll(".params td.name").length).toBe(2); // k and x
  });

  it("runs the model with the edited values and shows the result", async () => {
    const el = await mount();
    const k = el.querySelector("input[aria-label=k]") as HTMLInputElement;
    await act(async () => {
      k.value = "3";
      k.dispatchEvent(new Event("input", { bubbles: true }));
    });
    await click(button("Run"));
    await settle(300);
    const post = calls.find((c) => c.method === "POST" && c.path === "/runs")!;
    expect(post.body.model.parameters).toEqual({ k: 3 });
    expect(post.body.solver.endTime).toBe(4);
    expect(el.querySelector(".run-status.completed")).not.toBeNull();
    expect(el.querySelector("[data-testid=plot]")!.textContent).toBe("x,y");
    expect(el.textContent).toContain("1 event");
    expect(calls.some((c) => c.path.startsWith("/runs/r1/result?format=json&maxPoints="))).toBe(true);
  });

  it("sends the advanced solver settings and blocks values out of range", async () => {
    const el = await mount();
    const input = (label: string) => el.querySelector(`.advanced input[aria-label="${label}"]`) as HTMLInputElement;
    const type = async (label: string, value: string) => {
      await act(async () => {
        input(label).value = value;
        input(label).dispatchEvent(new Event("input", { bubbles: true }));
      });
    };
    await type("Error tolerance", "0.5");
    expect(input("Error tolerance").className).toContain("invalid");
    expect((button("Run") as HTMLButtonElement).disabled).toBe(true);
    await type("Error tolerance", "1e-7");
    await type("Longest step", "0.25");
    const mode = el.querySelector(".advanced select") as HTMLSelectElement;
    await act(async () => {
      mode.value = "smooth";
      mode.dispatchEvent(new Event("change", { bubbles: true }));
    });
    await type("Smoothing width", "1e-3");
    await click(button("Run"));
    await settle(300);
    const post = calls.find((c) => c.method === "POST" && c.path === "/runs")!;
    expect(post.body.solver).toEqual({ errorTol: 1e-7, segment: 0.25, endTime: 4 });
    expect(post.body.model.piecewise).toBe("smooth");
    expect(post.body.model.smoothEps).toBe(1e-3);
  });

  it("shows the result against a chosen variable and as a 3D trajectory", async () => {
    const el = await mount();
    await click(button("Run"));
    await settle(300);
    const choose = async (selector: string, value: string) => {
      const select = el.querySelector(selector) as HTMLSelectElement;
      await act(async () => {
        select.value = value;
        select.dispatchEvent(new Event("change", { bubbles: true }));
      });
    };
    await choose("select[aria-label='plot type']", "xy");
    expect(el.querySelector("[data-testid=plot]")).toBeNull();
    expect(el.querySelector("[data-testid=xy]")!.textContent).toBe("y against x"); // the first variable is on the axis
    await choose("select[aria-label='x axis']", "time");
    expect(el.querySelector("[data-testid=xy]")!.textContent).toBe("x,y against time");
    await click(el.querySelector(".legend input") as HTMLElement); // untick x
    expect(el.querySelector("[data-testid=xy]")!.textContent).toBe("y against time");
    await choose("select[aria-label='plot type']", "3d");
    expect(el.querySelector(".legend")).toBeNull();
    expect(el.querySelector("[data-testid=p3d]")!.textContent).toBe("x,y,y"); // two variables only: z repeats the last
    await choose("select[aria-label='z variable']", "x");
    expect(el.querySelector("[data-testid=p3d]")!.textContent).toBe("x,y,x");
    await choose("select[aria-label='plot type']", "time");
    expect(el.querySelector("[data-testid=plot]")).not.toBeNull();
  });

  it("does not run a model with errors and lists them", async () => {
    checkReply = { ok: false, diagnostics: [{ severity: "error", line: 2, column: 1, message: "broken", source: "parser" }], info: null };
    const el = await mount();
    expect(el.querySelector("ul.diagnostics")!.textContent).toContain("broken");
    expect((button("Run") as HTMLButtonElement).disabled).toBe(true);
    expect(el.querySelector(".status-line .error-text")!.textContent).toContain("1 error");
  });

  it("offers the models of the text and checks and runs the chosen one", async () => {
    checkReply = { ok: true, diagnostics: [], info: { ...info, models: ["Decay", "Other"], complete: true } };
    const el = await mount();
    const select = el.querySelector("select[aria-label=model]") as HTMLSelectElement;
    expect(select.value).toBe("Other"); // the last model is simulated unless one is chosen
    expect(calls.filter((c) => c.path === "/check").at(-1)!.body.className).toBeUndefined();
    await act(async () => {
      select.value = "Decay";
      select.dispatchEvent(new Event("change", { bubbles: true }));
    });
    await settle();
    expect(calls.filter((c) => c.path === "/check").at(-1)!.body.className).toBe("Decay");
    await click(button("Run"));
    await settle(300);
    expect(calls.find((c) => c.method === "POST" && c.path === "/runs")!.body.model.className).toBe("Decay");
    // a name that does not exist (any more) is dropped
    checkReply = { ok: false, diagnostics: [], info: { name: "", description: "", hybrid: false, parameters: [], variables: [], models: ["Other"], complete: false } };
    const editor = el.querySelector("[data-testid=editor]") as HTMLTextAreaElement;
    await act(async () => {
      editor.value = "model Other end Other;";
      editor.dispatchEvent(new Event("input", { bubbles: true }));
    });
    await settle();
    await settle(600);
    expect(calls.filter((c) => c.path === "/check").at(-1)!.body.className).toBeUndefined();
  });

  it("loads another example and clears the result", async () => {
    const el = await mount();
    await click(button("Run"));
    await settle(300);
    expect(el.querySelector("[data-testid=plot]")).not.toBeNull();
    const select = el.querySelector("select[aria-label=example]") as HTMLSelectElement;
    await act(async () => {
      select.value = "two";
      select.dispatchEvent(new Event("change", { bubbles: true }));
    });
    expect((el.querySelector("[data-testid=editor]") as HTMLTextAreaElement).value).toContain("model Two");
    expect((el.querySelector("input[aria-label='end time']") as HTMLInputElement).value).toBe("2");
    expect(el.querySelector("[data-testid=plot]")).toBeNull();
  });

  it("toggles one plot per variable and variable selection", async () => {
    const el = await mount();
    await click(button("Run"));
    await settle(300);
    const boxes = el.querySelectorAll("ul.legend input");
    await click(boxes[1]); // hide y
    expect(el.querySelector("[data-testid=plot]")!.textContent).toBe("x");
    await click(el.querySelector(".result-bar input[type=checkbox]"));
    expect(el.querySelector("[data-testid=plot]")!.textContent).toContain("stacked");
  });

  it("explains a server that cannot be reached", async () => {
    vi.spyOn(globalThis, "fetch").mockRejectedValue(new TypeError("down"));
    const el = await mount();
    expect(el.querySelector(".banner.bad")!.textContent).toContain("cannot reach the server");
    expect(el.querySelector(".conn.bad")).not.toBeNull();
  });

  it("asks for the token when the server wants one", async () => {
    vi.spyOn(globalThis, "fetch").mockResolvedValue(json({ error: { code: "unauthorized", message: "a valid bearer token is required" } }, 401));
    const el = await mount();
    expect(el.querySelector(".banner.bad")!.textContent).toContain("token");
  });

  it("saves new connection settings", async () => {
    const el = await mount();
    await click(button("Server"));
    const url = el.querySelector(".connection input[type=url]") as HTMLInputElement;
    await act(async () => {
      url.value = "http://other:9";
      url.dispatchEvent(new Event("input", { bubbles: true }));
    });
    await act(async () => {
      (el.querySelector(".connection") as HTMLFormElement).dispatchEvent(new Event("submit", { bubbles: true, cancelable: true }));
    });
    expect(JSON.parse(localStorage.getItem("gensas.settings")!).serverUrl).toBe("http://other:9");
  });

  it("warns about an API version it does not know", async () => {
    vi.spyOn(globalThis, "fetch").mockImplementation(async (input) => {
      const path = new URL(String(input)).pathname;
      return path === "/version" ? json({ ...version, api: 7 }) : json(path === "/examples" ? examples : checkReply);
    });
    const el = await mount();
    expect(el.querySelector(".banner.warn")!.textContent).toContain("API version 7");
  });

  it("cancels a running run", async () => {
    const encoder = new TextEncoder();
    let finish: () => void = () => {};
    const stream = new ReadableStream({
      start(controller) {
        controller.enqueue(encoder.encode(`event: state\ndata: ${JSON.stringify(run("running"))}\n\n`));
        finish = () => {
          controller.enqueue(encoder.encode(`event: done\ndata: ${JSON.stringify(run("cancelled", { error: "cancelled" }))}\n\n`));
          controller.close();
        };
      },
    });
    const el = await mount();
    const original = (globalThis.fetch as any).getMockImplementation();
    vi.spyOn(globalThis, "fetch").mockImplementation(async (input, init) => {
      if (new URL(String(input)).pathname === "/runs/r1/events") return new Response(stream);
      return original(input, init);
    });
    await click(button("Run"));
    await settle(100);
    expect(button("Cancel")).not.toBeNull();
    await click(button("Cancel"));
    await settle(50);
    expect(calls.some((c) => c.method === "DELETE" && c.path === "/runs/r1")).toBe(true);
    await act(async () => finish());
    await settle(100);
    expect(el.querySelector(".run-status.cancelled")).not.toBeNull();
    expect(button("Run")).not.toBeNull();
    expect(el.querySelector("[data-testid=plot]")).toBeNull();
  });

  it("shows why a run failed", async () => {
    runFinal = run("failed", { error: "the solver failed", diagnostics: [{ severity: "error", line: 3, column: 0, message: "bad model", source: "checker" }] });
    const el = await mount();
    await click(button("Run"));
    await settle(300);
    expect(el.querySelector(".run-status.failed")!.textContent).toContain("the solver failed");
    expect(el.textContent).toContain("bad model");
  });
});
