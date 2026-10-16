import { defaultKeymap, history, historyKeymap, indentWithTab } from "@codemirror/commands";
import { type Diagnostic as CmDiagnostic, lintGutter, setDiagnostics } from "@codemirror/lint";
import { EditorState } from "@codemirror/state";
import { EditorView, drawSelection, highlightActiveLine, keymap, lineNumbers } from "@codemirror/view";
import { useEffect, useRef } from "preact/hooks";
import { modelicaSupport } from "./modelica";
import type { Diagnostic } from "./types";

interface Props {
  value: string;
  onChange: (text: string) => void;
  diagnostics: Diagnostic[];
  /** set to a line number to scroll there (and select it); the object identity triggers the jump */
  jumpTo?: { line: number } | null;
}

/** Positions of a diagnostic in the document, clamped to it. */
export function toCodeMirror(doc: EditorState["doc"], d: Diagnostic): CmDiagnostic {
  const line = doc.line(Math.min(Math.max(d.line, 1), doc.lines));
  let from = line.from;
  let to = line.to;
  if (d.line >= 1 && d.column >= 1) {
    from = Math.min(line.from + d.column - 1, line.to);
    const word = /^[A-Za-z_][A-Za-z_0-9]*|^\S/.exec(line.text.slice(from - line.from));
    to = word ? Math.min(from + word[0].length, line.to) : from;
    if (to <= from) to = Math.min(from + 1, doc.length);
  } else if (d.line < 1) {
    from = 0;
    to = Math.min(doc.length, doc.line(1).to);
  }
  return { from, to, severity: d.severity === "error" ? "error" : "warning", message: d.message, source: d.source };
}

export function Editor({ value, onChange, diagnostics, jumpTo }: Props) {
  const host = useRef<HTMLDivElement>(null);
  const view = useRef<EditorView | null>(null);
  const onChangeRef = useRef(onChange);
  onChangeRef.current = onChange;

  useEffect(() => {
    const v = new EditorView({
      parent: host.current!,
      state: EditorState.create({
        doc: value,
        extensions: [
          lineNumbers(),
          lintGutter(),
          history(),
          drawSelection(),
          highlightActiveLine(),
          keymap.of([...defaultKeymap, ...historyKeymap, indentWithTab]),
          modelicaSupport,
          EditorView.updateListener.of((u) => {
            if (u.docChanged) onChangeRef.current(u.state.doc.toString());
          }),
        ],
      }),
    });
    view.current = v;
    return () => v.destroy();
  }, []);

  // text replaced from outside (an example was loaded)
  useEffect(() => {
    const v = view.current;
    if (v && v.state.doc.toString() !== value) {
      v.dispatch({ changes: { from: 0, to: v.state.doc.length, insert: value } });
    }
  }, [value]);

  useEffect(() => {
    const v = view.current;
    if (!v) return;
    v.dispatch(setDiagnostics(v.state, diagnostics.map((d) => toCodeMirror(v.state.doc, d))));
  }, [diagnostics, value]);

  useEffect(() => {
    const v = view.current;
    if (!v || !jumpTo) return;
    const line = v.state.doc.line(Math.min(Math.max(jumpTo.line, 1), v.state.doc.lines));
    v.dispatch({ selection: { anchor: line.from, head: line.to }, scrollIntoView: true });
    v.focus();
  }, [jumpTo]);

  return <div class="editor" ref={host} />;
}
