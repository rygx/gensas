import { EditorState } from "@codemirror/state";
import { describe, expect, it } from "vitest";
import { toCodeMirror } from "./Editor";

const doc = EditorState.create({ doc: "model M\n  Real x;\nequation\n  der(x) = -y;\nend M;\n" }).doc;

describe("toCodeMirror", () => {
  it("marks the identifier at the column", () => {
    const d = toCodeMirror(doc, { severity: "error", line: 4, column: 13, message: "m", source: "checker" });
    expect(doc.sliceString(d.from, d.to)).toBe("y");
    expect(d.severity).toBe("error");
  });

  it("marks the whole line when the column is unknown", () => {
    const d = toCodeMirror(doc, { severity: "warning", line: 2, column: 0, message: "m", source: "parser" });
    expect(doc.sliceString(d.from, d.to)).toBe("  Real x;");
    expect(d.severity).toBe("warning");
  });

  it("clamps a line past the end and handles line 0", () => {
    const far = toCodeMirror(doc, { severity: "error", line: 99, column: 1, message: "m", source: "x" });
    expect(far.from).toBeLessThanOrEqual(doc.length);
    const none = toCodeMirror(doc, { severity: "error", line: 0, column: 0, message: "m", source: "x" });
    expect(none.from).toBe(0);
  });
});
