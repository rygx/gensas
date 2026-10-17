export interface Settings {
  serverUrl: string;
  token: string;
}

const KEY = "gensas.settings";
const SOURCE_KEY = "gensas.source";

export function defaultServerUrl(): string {
  if (!import.meta.env.DEV && location.protocol.startsWith("http")) return location.origin; // served by the GenSAS server itself
  return "http://127.0.0.1:8764";
}

export function loadSettings(): Settings {
  const fallback = { serverUrl: defaultServerUrl(), token: "" };
  try {
    const raw = localStorage.getItem(KEY);
    return raw ? { ...fallback, ...JSON.parse(raw) } : fallback;
  } catch {
    return fallback;
  }
}

export function saveSettings(s: Settings): void {
  try {
    localStorage.setItem(KEY, JSON.stringify(s));
  } catch {
    /* storage may be unavailable */
  }
}

export function loadSource(): { source: string; endTime: number } | null {
  try {
    const raw = localStorage.getItem(SOURCE_KEY);
    const saved = raw ? JSON.parse(raw) : null;
    return saved && typeof saved.source === "string" && saved.source.trim() !== "" ? saved : null;
  } catch {
    return null;
  }
}

export function saveSource(source: string, endTime: number): void {
  if (source.trim() === "") return; // nothing worth keeping; the first example is shown instead
  try {
    localStorage.setItem(SOURCE_KEY, JSON.stringify({ source, endTime }));
  } catch {
    /* storage may be unavailable */
  }
}
