import { describe, expect, it } from "vitest";
import { DEFAULT_VIEW, clampView, cubeEdges, extent, formatTick, niceTicks, normalize, project, trailColor } from "./phase";

describe("extent", () => {
  it("ignores gaps and non-finite values", () => {
    expect(extent([[1, null, 3], [2, NaN, -1, Infinity]])).toEqual({ min: -1, max: 3 });
  });
  it("widens a constant and an empty series", () => {
    expect(extent([[2, 2]])).toEqual({ min: 1.8, max: 2.2 });
    expect(extent([[0, 0]])).toEqual({ min: -1, max: 1 });
    expect(extent([[]])).toEqual({ min: 0, max: 1 });
  });
});

describe("niceTicks", () => {
  it("uses steps of 1, 2 or 5 times a power of ten inside the range", () => {
    expect(niceTicks(0, 10, 5)).toEqual([0, 2, 4, 6, 8, 10]);
    expect(niceTicks(-2.3, 2.6, 5)).toEqual([-2, -1, 0, 1, 2]);
    expect(niceTicks(0, 0.0007, 7)).toEqual([0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007]);
  });
  it("copes with an empty range", () => {
    expect(niceTicks(3, 3)).toEqual([3]);
  });
});

describe("formatTick", () => {
  it("is short for normal and extreme numbers", () => {
    expect(formatTick(0)).toBe("0");
    expect(formatTick(2.5)).toBe("2.5");
    expect(formatTick(0.30000000000000004)).toBe("0.3");
    expect(formatTick(1e-5)).toBe("1.0e-5");
    expect(formatTick(250000)).toBe("2.5e5");
  });
});

describe("3D projection", () => {
  it("maps the range onto the cube", () => {
    expect(normalize(5, { min: 0, max: 10 })).toBe(0);
    expect(normalize(0, { min: 0, max: 10 })).toBe(-1);
    expect(normalize(10, { min: 0, max: 10 })).toBe(1);
  });

  it("looks along y at yaw 0 and pitch 0, with z up", () => {
    const view = { yaw: 0, pitch: 0, zoom: 1 };
    expect(project([0.5, 0.2, 0.8], view)).toEqual([0.5, 0.8, 0.2]);
  });

  it("turns about the vertical axis and tilts", () => {
    const [x, y, depth] = project([1, 0, 0], { yaw: Math.PI / 2, pitch: 0, zoom: 1 });
    expect(x).toBeCloseTo(0);
    expect(y).toBeCloseTo(0);
    expect(depth).toBeCloseTo(1);
    const tilted = project([0, 0, 1], { yaw: 0, pitch: Math.PI / 2, zoom: 2 });
    expect(tilted[0]).toBeCloseTo(0);
    expect(tilted[1]).toBeCloseTo(0); // looking down the z axis
    expect(tilted[2]).toBeCloseTo(2);
  });

  it("limits tilt and zoom", () => {
    expect(clampView({ yaw: 9, pitch: 5, zoom: 100 })).toEqual({ yaw: 9, pitch: Math.PI / 2, zoom: 6 });
    expect(clampView({ yaw: 0, pitch: -5, zoom: 0 })).toEqual({ yaw: 0, pitch: -Math.PI / 2, zoom: 0.3 });
    expect(clampView(DEFAULT_VIEW)).toEqual(DEFAULT_VIEW);
  });

  it("has the twelve edges of a cube", () => {
    expect(cubeEdges().length).toBe(12);
  });

  it("colors the start and end of a run differently", () => {
    expect(trailColor(0)).not.toBe(trailColor(1));
    expect(trailColor(-3)).toBe(trailColor(0));
  });
});
