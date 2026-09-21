# Codex Development Rules

These instructions apply to automated coding agents working in this repository. Read `DEVELOPMENT.md`, `README.md`, and the relevant files under `validation/` before making nontrivial changes.

## Priority

1. Follow the user's explicit task.
2. Preserve established physics and numerical semantics unless the task explicitly changes them.
3. Make the smallest coherent change that satisfies the task.
4. Do not silently broaden scope.

## Git safety

- Start by inspecting `git status`, the current branch, and the latest commit.
- Never discard, overwrite, reset, or rewrite user changes unless explicitly asked.
- Do not use destructive Git operations such as `git reset --hard`, `git clean -fd`, forced checkout, or force push without explicit permission.
- Do not commit directly to `main`. Use a short-lived branch named according to `DEVELOPMENT.md`.
- Do not push, merge, tag, publish a release, or change repository settings unless the user explicitly requests that action.
- Keep one scientific / numerical idea per branch and PR.

## Scope classification

Before editing, classify the task as one or more of:

- documentation;
- engineering / performance;
- numerical method;
- physics definition;
- analysis / observable;
- validation infrastructure.

State internally which class applies and choose tests accordingly.

A performance optimization is not automatically a physics change. A changed estimator, sampling distribution, trajectory termination rule, or physical boundary is.

## Current protected semantics

Unless explicitly requested, do not change:

- the incident ensemble / gravitational-focusing sampling law;
- the 1100 AU incident reference;
- the 2-Rsun injection checkpoint;
- the 1-Rsun numerical/analytic matching surface;
- the default 1100-Rsun bound-orbit removal prescription;
- the distinction between physical escape, outer-orbit removal, and evaporation;
- independent fixed-incident capture normalization;
- complete-history residence accounting;
- the exclusion of incomplete computational prefixes from production occupation;
- the radial-grid definition;
- the schema/provenance checks that prevent incompatible capture and transport products from being combined.

If a task requires changing one of these, treat it as a physics/numerical change: update tests, documentation, metadata/schema decisions, and validation together.

## Editing rules

- Prefer extending existing abstractions over creating parallel implementations.
- Preserve backward compatibility when it is cheap and unambiguous; reject ambiguous configuration rather than silently guessing.
- Do not silently reinterpret units.
- Numerical configuration that can affect results must be recorded in metadata.
- Runtime diagnostics must not be confused with physical uncertainty.
- Generated simulation outputs, build directories, raw validation runs, and local batch artifacts must not be committed.
- Do not change unrelated formatting or refactor unrelated code in the same task.

## Scientific correctness

- Distinguish deterministic regression tests from Monte Carlo convergence evidence.
- Do not claim convergence merely because two noisy estimates overlap.
- Do not classify a truncated / censored pilot as a production physics result.
- For nonlinear observables, preserve the existing uncertainty / bias treatment unless the task explicitly changes the estimator.
- Rare numerical failures in fixed-incident capture may be bounded explicitly; unresolved captured-transport failures must not be silently replaced by extra successful histories in production.
- Any weighted / importance-sampled transport requires a corresponding weighted occupation and uncertainty estimator; do not insert forced-collision histories into the current equal-weight estimator.

## Testing

During development, run the smallest focused test that exercises the change.

Before declaring a code-changing task complete, normally run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DCODE_COVERAGE=OFF
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

If an existing configured Release build is available and compatible, reuse it.

Additional expectations:

- MPI changes: run the relevant multi-rank regression.
- Python analysis changes: run `python3 scripts/test_analyze_point.py` and any directly affected script tests.
- numerical method changes: run a reference/convergence comparison in addition to unit tests.
- physics changes: run the corresponding protocol under `validation/`.
- documentation-only changes do not require an expensive rebuild unless they alter executable examples or workflow commands.

If a required test cannot be run, say exactly which test was not run and why. Never report a test as passed without executing it.

## Output / provenance

When adding a configuration knob that can alter numerical results:

1. validate its type, range, and disabled/enabled semantics;
2. make ambiguity an error rather than an implicit fallback;
3. record the effective value in metadata;
4. include the effective value in compatibility checks when products must match;
5. add a regression test for the configuration contract.

When preserving an old schema, only infer fields that are logically determined by the old metadata. Do not fabricate unavailable runtime diagnostics.

## Completion report

At the end of a task, report succinctly:

- branch and base commit;
- files changed;
- what behavior changed and what did not;
- exact tests run and their results;
- any remaining limitation;
- whether existing production results require regeneration.

Do not proceed to a second optimization or unrelated cleanup unless the user asks for it.
