/** Server-sent events: splits a text stream into {event, data} messages. */
export interface SseMessage {
  event: string;
  data: string;
}

export class SseParser {
  private buffer = "";

  /** Feeds a piece of the stream and returns the messages it completed. */
  push(chunk: string): SseMessage[] {
    this.buffer += chunk.replace(/\r\n/g, "\n");
    const out: SseMessage[] = [];
    let end: number;
    while ((end = this.buffer.indexOf("\n\n")) >= 0) {
      const block = this.buffer.slice(0, end);
      this.buffer = this.buffer.slice(end + 2);
      let event = "message";
      const data: string[] = [];
      for (const line of block.split("\n")) {
        if (line.startsWith(":")) continue; // comment (keep-alive)
        if (line.startsWith("event:")) event = line.slice(6).trim();
        else if (line.startsWith("data:")) data.push(line.slice(5).replace(/^ /, ""));
      }
      if (data.length > 0) out.push({ event, data: data.join("\n") });
    }
    return out;
  }
}
