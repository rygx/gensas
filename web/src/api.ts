import { SseParser } from "./sse";
import type { CheckResponse, Example, ModelSpec, ResultJson, Run, RunRequest, Version } from "./types";

/** The server answered with an error (`code` and `message` come from its error object). */
export class ApiError extends Error {
  constructor(
    public status: number,
    public code: string,
    message: string,
  ) {
    super(message);
  }
}

/** The server could not be reached. */
export class NetworkError extends Error {}

export interface ResultOptions {
  variables?: string[];
  maxPoints?: number;
}

export interface Columns {
  samples: number;
  /** time first, then the variables and discrete variables in the order of the run summary */
  columns: Float64Array[];
}

export class Api {
  constructor(
    public baseUrl: string,
    public token: string = "",
  ) {}

  private url(path: string): string {
    return this.baseUrl.replace(/\/+$/, "") + path;
  }

  private headers(json = false): Record<string, string> {
    const h: Record<string, string> = {};
    if (this.token) h["Authorization"] = "Bearer " + this.token;
    if (json) h["Content-Type"] = "application/json";
    return h;
  }

  private async send(method: string, path: string, body?: unknown, signal?: AbortSignal): Promise<Response> {
    let response: Response;
    try {
      response = await fetch(this.url(path), {
        method,
        headers: this.headers(body !== undefined),
        body: body === undefined ? undefined : JSON.stringify(body),
        signal,
      });
    } catch (err) {
      if ((err as Error).name === "AbortError") throw err;
      throw new NetworkError(`cannot reach the server at ${this.baseUrl}`);
    }
    if (!response.ok) {
      let code = "error";
      let message = response.statusText || "request failed";
      try {
        const e = (await response.json()).error;
        code = e.code;
        message = e.message;
      } catch {
        /* not an API error body */
      }
      throw new ApiError(response.status, code, message);
    }
    return response;
  }

  private async json<T>(method: string, path: string, body?: unknown): Promise<T> {
    return (await this.send(method, path, body)).json() as Promise<T>;
  }

  version(): Promise<Version> {
    return this.json("GET", "/version");
  }
  examples(): Promise<Example[]> {
    return this.json("GET", "/examples");
  }
  check(model: ModelSpec): Promise<CheckResponse> {
    return this.json("POST", "/check", model);
  }
  createRun(request: Omit<RunRequest, "kind">): Promise<Run> {
    return this.json("POST", "/runs", { kind: "modelica", ...request });
  }
  getRun(id: string): Promise<Run> {
    return this.json("GET", "/runs/" + encodeURIComponent(id));
  }
  listRuns(): Promise<Run[]> {
    return this.json("GET", "/runs");
  }
  /** Stops a waiting or running run. Returns the run, or null if it was finished and has been removed. */
  async cancelRun(id: string): Promise<Run | null> {
    const r = await this.send("DELETE", "/runs/" + encodeURIComponent(id));
    return r.status === 204 ? null : ((await r.json()) as Run);
  }

  private resultQuery(format: string, o: ResultOptions): string {
    const q = new URLSearchParams({ format });
    if (o.variables && o.variables.length > 0) q.set("variables", o.variables.join(","));
    if (o.maxPoints) q.set("maxPoints", String(o.maxPoints));
    return "/result?" + q.toString();
  }

  async result(id: string, o: ResultOptions = {}): Promise<ResultJson> {
    const r = await this.send("GET", "/runs/" + encodeURIComponent(id) + this.resultQuery("json", o));
    return (await r.json()) as ResultJson;
  }

  /** The samples as raw columns (much smaller than JSON for big results). `names` must be the columns after time. */
  async resultColumns(id: string, o: ResultOptions = {}): Promise<Columns> {
    const r = await this.send("GET", "/runs/" + encodeURIComponent(id) + this.resultQuery("f64", o));
    return decodeColumns(await r.arrayBuffer(), Number(r.headers.get("X-Gensas-Samples")));
  }

  /** Downloads a result file (csv or mat) as a Blob. */
  async resultFile(id: string, format: "csv" | "mat"): Promise<Blob> {
    return (await this.send("GET", "/runs/" + encodeURIComponent(id) + this.resultQuery(format, {}))).blob();
  }

  /**
   * Follows a run until it ends, calling `onUpdate` with every state and progress change. Resolves with the
   * finished run. Uses the event stream; falls back to polling if the stream cannot be read.
   */
  async follow(id: string, onUpdate: (run: Run) => void, signal?: AbortSignal): Promise<Run> {
    try {
      return await this.followStream(id, onUpdate, signal);
    } catch (err) {
      if ((err as Error).name === "AbortError") throw err;
      return this.followPolling(id, onUpdate, signal);
    }
  }

  private async followStream(id: string, onUpdate: (run: Run) => void, signal?: AbortSignal): Promise<Run> {
    const response = await this.send("GET", "/runs/" + encodeURIComponent(id) + "/events", undefined, signal);
    if (!response.body) throw new Error("no stream");
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    const parser = new SseParser();
    let last = null as Run | null;
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      for (const m of parser.push(decoder.decode(value, { stream: true }))) {
        if (m.event === "progress" && last) {
          last = { ...last, progress: JSON.parse(m.data) };
          onUpdate(last);
        } else if (m.event === "state" || m.event === "done") {
          last = JSON.parse(m.data) as Run;
          onUpdate(last);
          if (m.event === "done") return last;
        }
      }
    }
    throw new Error("the event stream ended early");
  }

  private async followPolling(id: string, onUpdate: (run: Run) => void, signal?: AbortSignal): Promise<Run> {
    for (;;) {
      if (signal?.aborted) throw new DOMException("aborted", "AbortError");
      const run = await this.getRun(id);
      onUpdate(run);
      if (run.state !== "queued" && run.state !== "running") return run;
      await new Promise((resolve) => setTimeout(resolve, 300));
    }
  }
}

export function decodeColumns(buffer: ArrayBuffer, samples: number): Columns {
  const all = new Float64Array(buffer);
  const columns: Float64Array[] = [];
  for (let start = 0; start + samples <= all.length && samples > 0; start += samples) {
    columns.push(all.subarray(start, start + samples));
  }
  return { samples, columns };
}
