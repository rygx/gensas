import { describe, expect, it } from "vitest";
import { advancedErrors, advancedRequest, isAdvancedChanged, piecewiseOf } from "./advanced";

describe("advanced settings", () => {
  it("sends nothing by default", () => {
    expect(advancedRequest({})).toEqual({ solver: {}, model: {} });
    expect(isAdvancedChanged({})).toBe(false);
    expect(piecewiseOf({})).toBe("events");
  });

  it("turns the typed values into request fields", () => {
    const r = advancedRequest({ segment: "0.5", level: "20", alphaTol: "1e-4", errorTol: "1e-7", outputStep: "0.02", piecewise: "smooth", smoothEps: "1e-3" });
    expect(r.solver).toEqual({ segment: 0.5, level: 20, alphaTol: 1e-4, errorTol: 1e-7 });
    expect(r.model).toEqual({ piecewise: "smooth", smoothEps: 1e-3 });
    expect(r.outputStep).toBe(0.02);
    expect(isAdvancedChanged({ level: "20" })).toBe(true);
  });

  it("flags values outside the engine's ranges and does not send them", () => {
    const edits = { segment: "100", level: "2.5", errorTol: "0.5", outputStep: "x", alphaTol: "20" };
    const e = advancedErrors(edits);
    expect(Object.keys(e).sort()).toEqual(["alphaTol", "errorTol", "level", "outputStep", "segment"]);
    expect(advancedRequest(edits)).toEqual({ solver: {}, model: {} });
  });

  it("limits the step shrink tolerance by the longest step", () => {
    expect(advancedErrors({ segment: "0.1", alphaTol: "0.5" }).alphaTol).toContain("0.1");
    expect(advancedErrors({ segment: "0.1", alphaTol: "0.05" }).alphaTol).toBeUndefined();
  });

  it("ignores the smoothing width while the exact mode is chosen", () => {
    expect(advancedErrors({ smoothEps: "-1" })).toEqual({});
    expect(advancedErrors({ piecewise: "smooth", smoothEps: "-1" }).smoothEps).toBe("must be positive");
    expect(advancedRequest({ smoothEps: "1e-3" }).model).toEqual({});
  });
});
