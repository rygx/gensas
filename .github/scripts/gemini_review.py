#!/usr/bin/env python3
"""
Gemini AI code reviewer for pull requests (GenSAS).

Called by .github/workflows/gemini-review.yml on pull_request events.
Fetches the diff for reviewable files only, sends it to Gemini, and posts the
result as a PR comment (replacing any previous Gemini review).

Environment variables (injected by the workflow):
  GITHUB_TOKEN       - for posting comments
  GEMINI_API_KEY     - Google AI Studio key (free tier supported)
  GITHUB_REPOSITORY  - owner/repo
  PR_NUMBER          - pull request number
  BASE_SHA           - base commit SHA of the PR
  HEAD_SHA           - head commit SHA of the PR

The review is advisory: if every model is out of quota, a soft notice is posted
and the job exits 0 so it never blocks a merge.
"""

import os
import re
import subprocess
import sys
import time
import traceback

import google.generativeai as genai
from google.api_core.exceptions import ResourceExhausted
from github import Github

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Models are tried in order, best first; a model that is out of quota or fails is skipped and
# the next one is used, and every model is attempted once per round before any backoff.
# Verify IDs against https://ai.google.dev/gemini-api/docs/models when changing this list.
MODELS = [
    "gemini-3.8-flash",  # newest stable Flash (checked on ai.google.dev, Oct 2026)
    "gemini-3.6-flash",
    "gemini-3.5-flash-lite",
    "gemini-3.1-flash-lite",
    "gemini-2.5-flash",
    "gemini-2.5-flash-lite",
    "gemma-4-26b-a4b-it",
    "gemma-4-31b-it",
]

# Hard cap on diff characters sent to Gemini (free-tier per-request token limits).
MAX_DIFF_CHARS = 60_000

# One "round" tries every model once. Backoff between rounds: 30s, 60s, 120s.
MAX_ROUNDS = 4
INITIAL_ROUND_BACKOFF_SECONDS = 30

# Extensions worth reviewing in this repository (C++ / Bazel / CI / docs).
REVIEW_EXTENSIONS = {
    ".cpp", ".cc", ".cxx", ".c", ".h", ".hpp", ".hxx",  # C/C++ sources and headers
    ".bazel", ".bzl",                           # BUILD files and Starlark
    ".yml", ".yaml",                            # GitHub Actions
    ".sh", ".py",                               # scripts
    ".mo",                                      # Modelica models
    ".md",                                      # documentation
    ".txt",                                     # token/opinfo tables (X-macro lists)
}

# Exact filenames to skip even when their extension is allowed.
SKIP_FILES = {"MODULE.bazel.lock", "LICENSE.md"}

REVIEW_HEADER = "## Gemini Code Review 🤖"

_DIFF_HEADER_RE = re.compile(r"^diff --git a/.+ b/(.+)$")

PROMPT_TEMPLATE = """\
You are an expert code reviewer for GenSAS, a C++17 generic simulation tool based on
semi-analytical solutions (SAS / Pade approximation). It parses a subset of Modelica,
solves ODE/DAE models segment by segment, and also runs AC power flow on PSAT cases.
It is built with Bazel (WORKSPACE, rules_foreign_cc) on Linux, macOS, and Windows
(MSYS2 MinGW-w64), uses Armadillo/SuperLU/OpenBLAS, and keeps all three platforms working.

Review the following git diff and provide concise, actionable feedback.

Focus on:
- Bugs, logic errors, undefined behavior, memory/ownership problems (raw pointers, leaks)
- Numerical correctness and robustness (tolerances, division by zero, NaN handling)
- Parser/lexer correctness and quality of diagnostics
- Cross-platform issues (Linux, macOS, Windows/MinGW): paths, line endings, 32/64-bit types
- Bazel/BUILD correctness (deps, visibility, select() on platforms)
- Test quality: do tests actually pin the intended behavior?
- Naming, readability, and maintainability

For documentation files (.md): technical accuracy against the code, completeness, clarity.

Format your response exactly as follows:

## Gemini Code Review 🤖

### Summary
[1-2 sentences]

### Issues
[Each issue: `file:line - 🔴 Critical / 🟠 Major / 🟡 Minor - description and fix`]
[Write "None found." if there are no issues]

### Suggestions
[Optional improvements that are not bugs. Write "None." if nothing to add]

### Looks good ✅
[What is done well]

Be concise. No padding or filler.

---
{diff}
"""


# ---------------------------------------------------------------------------
# Diff helpers
# ---------------------------------------------------------------------------


def _should_review(filename: str) -> bool:
    ext = os.path.splitext(filename)[1].lower()
    basename = os.path.basename(filename)
    return (ext in REVIEW_EXTENSIONS or basename in {"BUILD", "WORKSPACE"}) and basename not in SKIP_FILES


def get_filtered_diff(base_sha: str, head_sha: str) -> tuple[str, bool]:
    """Return (diff, truncated) for reviewable files only."""
    names_result = subprocess.run(
        ["git", "diff", "--name-only", f"{base_sha}...{head_sha}"],
        capture_output=True,
        text=True,
        check=True,
    )
    all_files = [f.strip() for f in names_result.stdout.strip().splitlines() if f.strip()]
    reviewable = [f for f in all_files if _should_review(f)]

    if not reviewable:
        return "", False

    diff_result = subprocess.run(
        ["git", "diff", f"{base_sha}...{head_sha}", "--"] + reviewable,
        capture_output=True,
        text=True,
        check=True,
    )
    diff = diff_result.stdout

    if len(diff) > MAX_DIFF_CHARS:
        return diff[:MAX_DIFF_CHARS], True
    return diff, False


# ---------------------------------------------------------------------------
# Gemini
# ---------------------------------------------------------------------------


def _call_model_once(model_name: str, diff: str, truncated: bool) -> str:
    """Single attempt at one model; raises on failure."""
    model = genai.GenerativeModel(model_name)

    notice = (
        f"\n\n> ⚠️ Diff was truncated to {MAX_DIFF_CHARS:,} chars. "
        "Some files may not have been reviewed.\n"
        if truncated
        else ""
    )

    prompt = PROMPT_TEMPLATE.format(diff=diff) + notice
    response = model.generate_content(prompt)
    try:
        return response.text
    except ValueError as exc:
        # Empty candidates, usually a safety-filter block: fall through to the next model.
        raise RuntimeError(
            f"Response blocked or empty (prompt_feedback={response.prompt_feedback}): {exc}"
        ) from exc


def call_gemini(diff: str, truncated: bool) -> tuple[str, str]:
    """Round-based rotation over MODELS; returns (review_text, used_model_name)."""
    round_backoff = INITIAL_ROUND_BACKOFF_SECONDS

    for round_num in range(1, MAX_ROUNDS + 1):
        for model_name in MODELS:
            print(f"[Round {round_num}] Trying {model_name}...", file=sys.stderr)
            try:
                review = _call_model_once(model_name, diff, truncated)
                return review, model_name
            except ResourceExhausted:
                print(f"[{model_name}] Quota exhausted - trying next model.", file=sys.stderr)
            except Exception as exc:
                print(
                    f"[{model_name}] API error ({exc.__class__.__name__}: {exc}) - trying next model.",
                    file=sys.stderr,
                )

        if round_num < MAX_ROUNDS:
            print(
                f"All models failed in round {round_num}/{MAX_ROUNDS} - "
                f"backing off {round_backoff}s...",
                file=sys.stderr,
            )
            time.sleep(round_backoff)
            round_backoff *= 2

    raise RuntimeError(f"All {len(MODELS)} models exhausted after {MAX_ROUNDS} rounds.")


# ---------------------------------------------------------------------------
# GitHub comment
# ---------------------------------------------------------------------------


def post_or_update_comment(review: str) -> None:
    gh = Github(os.environ["GITHUB_TOKEN"])
    repo = gh.get_repo(os.environ["GITHUB_REPOSITORY"])
    pr = repo.get_pull(int(os.environ["PR_NUMBER"]))

    # Replace any previous Gemini review comment to keep the PR tidy.
    for comment in pr.get_issue_comments():
        if comment.body.startswith(REVIEW_HEADER):
            comment.delete()

    pr.create_issue_comment(review)
    print(f"Review posted on PR #{os.environ['PR_NUMBER']}.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    if not os.environ.get("GEMINI_API_KEY"):
        print("GEMINI_API_KEY is not set; skipping review (advisory job).")
        return

    genai.configure(api_key=os.environ["GEMINI_API_KEY"])

    base_sha = os.environ["BASE_SHA"]
    head_sha = os.environ["HEAD_SHA"]

    print(f"Diffing {base_sha[:8]}...{head_sha[:8]}")
    diff, truncated = get_filtered_diff(base_sha, head_sha)

    if not diff.strip():
        print("No reviewable changes in this diff (all changed files are excluded).")
        return

    print(
        f"Sending {len(diff):,} chars to Gemini "
        f"({len(MODELS)} models, up to {MAX_ROUNDS} rounds)..."
    )

    try:
        review_text, used_model = call_gemini(diff, truncated)
    except RuntimeError:
        gh = Github(os.environ["GITHUB_TOKEN"])
        repo = gh.get_repo(os.environ["GITHUB_REPOSITORY"])
        pr = repo.get_pull(int(os.environ["PR_NUMBER"]))
        pr.create_issue_comment(
            f"{REVIEW_HEADER}\n\n"
            "> ⏸️ Gemini review skipped - free-tier quota exhausted on all models "
            f"({', '.join(MODELS)}) across {MAX_ROUNDS} rounds. "
            "No action required; this does not block merge."
        )
        print("All models/rounds exhausted - soft notice posted, job exits 0.")
        return

    review_with_attribution = review_text.rstrip() + f"\n\n---\n*Review by `{used_model}`*"
    post_or_update_comment(review_with_attribution)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        traceback.print_exc()
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
