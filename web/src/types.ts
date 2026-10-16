import type { components } from "./api-types";

type S = components["schemas"];
export type Version = S["Version"];
export type Example = S["Example"];
export type ModelSpec = S["ModelSpec"];
export type SolverSpec = S["SolverSpec"];
export type RunRequest = S["RunRequest"];
export type Diagnostic = S["Diagnostic"];
export type ModelInfo = S["ModelInfo"];
export type ModelVariable = S["ModelVariable"];
export type ModelParameter = S["ModelParameter"];
export type CheckResponse = S["CheckResponse"];
export type Run = S["Run"];
export type ResultJson = S["ResultJson"];

export type RunState = Run["state"];
export const FINISHED: RunState[] = ["completed", "terminated", "failed", "cancelled"];
export const HAS_RESULT: RunState[] = ["completed", "terminated"];
