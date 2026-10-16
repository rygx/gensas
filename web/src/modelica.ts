import { HighlightStyle, StreamLanguage, syntaxHighlighting } from "@codemirror/language";
import { tags } from "@lezer/highlight";

const KEYWORDS = new Set(
  ("model end equation algorithm parameter constant input output discrete when elsewhen then if elseif else for in loop while " +
    "function return annotation each and or not true false within import extends class block record package type connector").split(" "),
);
const TYPES = new Set(["Real", "Integer", "Boolean", "String"]);
const BUILTINS = new Set(
  ("der pre reinit terminate assert initial edge change noEvent abs min max sign sqrt exp log log10 sin cos tan asin acos atan atan2 " +
    "sinh cosh tanh floor ceil integer mod rem time").split(" "),
);

/** A small Modelica tokenizer for highlighting; the server is the authority on what is valid. */
const modelica = StreamLanguage.define<{ comment: boolean }>({
  startState: () => ({ comment: false }),
  token(stream, state) {
    if (state.comment) {
      if (stream.match(/^.*?\*\//)) state.comment = false;
      else stream.skipToEnd();
      return "comment";
    }
    if (stream.eatSpace()) return null;
    if (stream.match("//")) {
      stream.skipToEnd();
      return "comment";
    }
    if (stream.match("/*")) {
      state.comment = true;
      return "comment";
    }
    if (stream.match(/^"(?:[^"\\]|\\.)*"?/)) return "string";
    if (stream.match(/^\d+\.?\d*(?:[eE][+-]?\d+)?|^\.\d+/)) return "number";
    if (stream.match(/^[A-Za-z_][A-Za-z_0-9]*/)) {
      const word = stream.current();
      if (KEYWORDS.has(word)) return "keyword";
      if (TYPES.has(word)) return "typeName";
      if (BUILTINS.has(word)) return "function";
      return "variableName";
    }
    stream.next();
    return null;
  },
});

const style = HighlightStyle.define([
  { tag: tags.keyword, color: "var(--syn-keyword)", fontWeight: "600" },
  { tag: tags.typeName, color: "var(--syn-type)" },
  { tag: tags.comment, color: "var(--syn-comment)", fontStyle: "italic" },
  { tag: tags.string, color: "var(--syn-string)" },
  { tag: tags.number, color: "var(--syn-number)" },
  { tag: tags.function(tags.variableName), color: "var(--syn-function)" },
]);

export const modelicaSupport = [modelica, syntaxHighlighting(style)];
