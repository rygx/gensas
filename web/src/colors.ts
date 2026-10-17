export const PALETTE = ["#2f6fdb", "#d9480f", "#2b8a3e", "#9c36b5", "#c2255c", "#0b7285", "#e67700", "#5f3dc4", "#495057", "#66a80f"];

export function colorOf(index: number): string {
  return PALETTE[index % PALETTE.length];
}
