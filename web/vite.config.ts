import { defineConfig } from "vitest/config";
import preact from "@preact/preset-vite";

export default defineConfig({
  plugins: [preact()],
  base: "./", // relative paths: the build can be hosted at any URL, or served by the GenSAS server
  test: { environment: "jsdom", include: ["src/**/*.test.ts"] },
});
