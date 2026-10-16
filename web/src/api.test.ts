import { afterEach, describe, expect, it, vi } from "vitest";
import { Api, ApiError, NetworkError, decodeColumns } from "./api";
import { SseParser } from "./sse";

function reply(body: unknown, status = 200, headers: Record<string, string> = {}): Response {
  const text = typeof body === "string" || body instanceof ArrayBuffer ? body : JSON.stringify(body);
  return new Response(status === 204 ? null : text, { status, headers });
}

afterEach(() => vi.restoreAllMocks());

describe("Api", () => {
  it("sends the token and the JSON body", async () => {
    const fetchMock = vi.spyOn(globalThis, "fetch").mockResolvedValue(reply({ ok: true, diagnostics: [], info: null }));
    const api = new Api("http://host:1/", "secret");
    await api.check({ source: "model M end M;" });
    const [url, init] = fetchMock.mock.calls[0];
    expect(url).toBe("http://host:1/check");
    expect((init!.headers as Record<string, string>).Authorization).toBe("Bearer secret");
    expect(JSON.parse(init!.body as string)).toEqual({ source: "model M end M;" });
  });

  it("turns error bodies into ApiError", async () => {
    vi.spyOn(globalThis, "fetch").mockResolvedValue(reply({ error: { code: "too_large", message: "big" } }, 413));
    await expect(new Api("http://h").check({ source: "x" })).rejects.toMatchObject({ status: 413, code: "too_large", message: "big" });
    await expect(new Api("http://h").check({ source: "x" })).rejects.toBeInstanceOf(ApiError);
  });

  it("reports an unreachable server", async () => {
    vi.spyOn(globalThis, "fetch").mockRejectedValue(new TypeError("failed"));
    await expect(new Api("http://nowhere").version()).rejects.toBeInstanceOf(NetworkError);
  });

  it("builds the result query", async () => {
    const fetchMock = vi.spyOn(globalThis, "fetch").mockResolvedValue(reply({}));
    await new Api("http://h").result("abc", { variables: ["x", "y"], maxPoints: 500 });
    expect(fetchMock.mock.calls[0][0]).toBe("http://h/runs/abc/result?format=json&variables=x%2Cy&maxPoints=500");
  });

  it("cancelRun returns null when the finished run was removed", async () => {
    vi.spyOn(globalThis, "fetch").mockResolvedValue(reply("", 204));
    expect(await new Api("http://h").cancelRun("r")).toBeNull();
  });

  it("follows the event stream to the done event", async () => {
    const run = (state: string) => ({ id: "r", state, progress: { time: 0, end: 1 } });
    const text =
      `event: state\ndata: ${JSON.stringify(run("running"))}\n\n` +
      `event: progress\ndata: ${JSON.stringify({ time: 0.5, end: 1 })}\n\n` +
      `event: done\ndata: ${JSON.stringify(run("completed"))}\n\n`;
    vi.spyOn(globalThis, "fetch").mockResolvedValue(reply(text));
    const seen: [string, number][] = [];
    const final = await new Api("http://h").follow("r", (r) => seen.push([r.state, r.progress.time]));
    expect(final.state).toBe("completed");
    expect(seen).toEqual([
      ["running", 0],
      ["running", 0.5],
      ["completed", 0],
    ]);
  });

  it("falls back to polling when the stream fails", async () => {
    let calls = 0;
    vi.spyOn(globalThis, "fetch").mockImplementation(async (url) => {
      if (String(url).endsWith("/events")) return reply({ error: { code: "x", message: "no" } }, 500);
      calls++;
      return reply({ id: "r", state: calls < 2 ? "running" : "completed", progress: { time: 0, end: 1 } });
    });
    const final = await new Api("http://h").follow("r", () => {});
    expect(final.state).toBe("completed");
    expect(calls).toBe(2);
  });
});

describe("decodeColumns", () => {
  it("splits the matrix into columns", () => {
    const data = new Float64Array([0, 1, 2, 10, 11, 12]);
    const { columns, samples } = decodeColumns(data.buffer, 3);
    expect(samples).toBe(3);
    expect(columns.map((c) => Array.from(c))).toEqual([
      [0, 1, 2],
      [10, 11, 12],
    ]);
  });
});

describe("SseParser", () => {
  it("handles messages split across chunks, comments and CRLF", () => {
    const p = new SseParser();
    expect(p.push("event: a\nda")).toEqual([]);
    expect(p.push("ta: 1\n\n: keep-alive\n\nevent: b\r\ndata: 2\r\n\r\n")).toEqual([
      { event: "a", data: "1" },
      { event: "b", data: "2" },
    ]);
  });
});
