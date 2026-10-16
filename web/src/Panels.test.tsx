import { act } from "preact/test-utils";
import { render } from "preact";
import { afterEach, describe, expect, it, vi } from "vitest";
import { DiagnosticList, type Edits, ParamsPanel, RunStatus, overrides, parseNumber } from "./Panels";
import type { ModelInfo, Run } from "./types";

const info: ModelInfo = {
  name: "M",
  description: "",
  hybrid: false,
  parameters: [
    { name: "k", type: "Real", unit: "1/s", description: "rate", value: 2, constant: false, line: 2 },
    { name: "c", type: "Real", unit: "", description: "", value: 3, constant: true, line: 3 },
  ],
  variables: [
    { name: "x", kind: "state", type: "Real", unit: "m", description: "", start: 1, input: false, line: 4 },
    { name: "y", kind: "algebraic", type: "Real", unit: "", description: "", start: null, input: false, line: 5 },
  ],
};

let host: HTMLElement | null = null;
function mount(vnode: preact.ComponentChild): HTMLElement {
  host = document.createElement("div");
  document.body.append(host);
  act(() => render(vnode, host!));
  return host;
}
afterEach(() => {
  if (host) act(() => render(null, host!));
  host?.remove();
  host = null;
});

describe("number fields", () => {
  it("parse numbers and refuse the rest", () => {
    expect(parseNumber("2.5")).toBe(2.5);
    expect(parseNumber("1e-3")).toBe(0.001);
    expect(parseNumber("")).toBeNull();
    expect(parseNumber("1e")).toBeNull();
    expect(parseNumber("abc")).toBeNull();
  });

  it("send only edits that are numbers", () => {
    expect(overrides({ a: "1", b: "", c: "x", d: "-2e1" })).toEqual({ a: 1, d: -20 });
  });
});

describe("ParamsPanel", () => {
  it("lists parameters and the states, but not algebraic variables", () => {
    const el = mount(<ParamsPanel info={info} params={{}} start={{}} setParams={() => {}} setStart={() => {}} />);
    expect(Array.from(el.querySelectorAll("td.name")).map((n) => n.textContent)).toEqual(["k", "c", "x"]);
    expect((el.querySelector("input[aria-label=c]") as HTMLInputElement).disabled).toBe(true);
    expect((el.querySelector("input[aria-label=k]") as HTMLInputElement).value).toBe("2");
  });

  it("reports an edit and offers to reset it", () => {
    const setParams = vi.fn();
    let params: Edits = {};
    const view = () => <ParamsPanel info={info} params={params} start={{}} setParams={setParams} setStart={() => {}} />;
    const el = mount(view());
    const input = el.querySelector("input[aria-label=k]") as HTMLInputElement;
    act(() => {
      input.value = "5";
      input.dispatchEvent(new Event("input", { bubbles: true }));
    });
    expect(setParams).toHaveBeenCalledWith({ k: "5" });
    params = { k: "5" };
    act(() => render(view(), el));
    expect(el.querySelector("tr.changed")).not.toBeNull();
    act(() => (el.querySelector("td.reset button") as HTMLButtonElement).click());
    expect(setParams).toHaveBeenLastCalledWith({});
  });

  it("marks a field that is not a number", () => {
    const el = mount(<ParamsPanel info={info} params={{ k: "1e" }} start={{}} setParams={() => {}} setStart={() => {}} />);
    expect(el.querySelector("input[aria-label=k]")!.className).toContain("invalid");
  });

  it("explains itself without model info", () => {
    const el = mount(<ParamsPanel info={null} params={{}} start={{}} setParams={() => {}} setStart={() => {}} />);
    expect(el.textContent).toContain("no errors");
  });
});

describe("DiagnosticList", () => {
  it("jumps to the line of a diagnostic", () => {
    const jump = vi.fn();
    const el = mount(
      <DiagnosticList
        items={[
          { severity: "error", line: 4, column: 3, message: "bad", source: "checker" },
          { severity: "warning", line: 0, column: 0, message: "general", source: "job" },
        ]}
        onJump={jump}
      />,
    );
    const buttons = el.querySelectorAll("button");
    act(() => (buttons[0] as HTMLButtonElement).click());
    expect(jump).toHaveBeenCalledWith(4);
    expect((buttons[1] as HTMLButtonElement).disabled).toBe(true);
    expect(el.querySelector("li.error")!.textContent).toContain("bad");
  });
});

describe("RunStatus", () => {
  const base: Run = { id: "r", kind: "modelica", state: "running", createdAt: 0, startedAt: 0, finishedAt: null, progress: { time: 1, end: 4 }, error: null, diagnostics: [], warnings: [], summary: null };
  it("shows progress while running", () => {
    const el = mount(<RunStatus run={base} />);
    expect((el.querySelector("progress") as HTMLProgressElement).value).toBe(0.25);
  });
  it("shows the error and warnings of a failed run", () => {
    const el = mount(<RunStatus run={{ ...base, state: "failed", error: "boom", warnings: ["clamped"] }} />);
    expect(el.textContent).toContain("failed");
    expect(el.textContent).toContain("boom");
    expect(el.textContent).toContain("clamped");
    expect(el.querySelector("progress")).toBeNull();
  });
  it("shows the message of terminate()", () => {
    const el = mount(<RunStatus run={{ ...base, state: "terminated", termination: { time: 0.5, message: "done" } }} />);
    expect(el.textContent).toContain("done");
  });
  it("renders nothing without a run", () => {
    const el = mount(<RunStatus run={null} />);
    expect(el.textContent).toBe("");
  });
});
