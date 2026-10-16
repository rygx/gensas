import { ensureSyntaxTree } from "@codemirror/language";
import { EditorState } from "@codemirror/state";
import { describe, expect, it } from "vitest";
import { modelicaSupport } from "./modelica";

function tokens(text: string): Map<string, string[]> {
  const state = EditorState.create({ doc: text, extensions: [modelicaSupport] });
  const tree = ensureSyntaxTree(state, state.doc.length, 5000)!;
  const found = new Map<string, string[]>();
  tree.iterate({
    enter(node) {
      if (node.name && node.from < node.to) {
        const list = found.get(node.name) ?? [];
        list.push(text.slice(node.from, node.to));
        found.set(node.name, list);
      }
    },
  });
  return found;
}

describe("Modelica language support", () => {
  it("recognises keywords, types, numbers, strings and both kinds of comments", () => {
    const t = tokens('model M "d"\n  /* block\n comment */ Real x(start = 1.5e-3); // line\nequation\n  der(x) = -x;\nend M;\n');
    const all = (name: string) => (t.get(name) ?? []).join("|");
    expect(all("keyword")).toContain("model");
    expect(all("keyword")).toContain("equation");
    expect(all("typeName")).toContain("Real");
    expect(all("number")).toContain("1.5e-3");
    expect(all("string")).toContain('"d"');
    expect(all("comment")).toContain("// line");
    expect(all("comment")).toContain("/* block");
    expect(all("comment")).toContain("comment */");
  });
});
