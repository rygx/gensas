import { beforeEach, describe, expect, it } from "vitest";
import { loadSettings, loadSource, saveSettings, saveSource } from "./settings";

beforeEach(() => localStorage.clear());

describe("settings", () => {
  it("round-trips the connection", () => {
    saveSettings({ serverUrl: "http://example:1", token: "t" });
    expect(loadSettings()).toEqual({ serverUrl: "http://example:1", token: "t" });
  });

  it("falls back to defaults on bad stored data", () => {
    localStorage.setItem("gensas.settings", "{not json");
    expect(loadSettings().token).toBe("");
  });

  it("keeps the model but never an empty one", () => {
    expect(loadSource()).toBeNull();
    saveSource("model M end M;", 3);
    expect(loadSource()).toEqual({ source: "model M end M;", endTime: 3 });
    saveSource("   ", 4);
    expect(loadSource()!.endTime).toBe(3);
  });
});
