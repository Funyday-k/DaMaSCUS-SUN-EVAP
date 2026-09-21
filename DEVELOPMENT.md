# Development and Production Policy

This repository is research software. The development process must preserve not only software correctness but also the traceability of every physics result to a specific numerical definition, source revision, configuration, and validation record.

## 1. Repository roles

- `main` is the stable integration branch. It should build, pass the automated test gate, and have no known unresolved physics or numerical regression.
- Development is performed on short-lived branches and merged through pull requests.
- Production campaigns must use an exact clean commit or, preferably, a release tag. Do not run paper production from an unspecified moving `main`.
- Large generated simulation outputs stay outside Git. Small manifests, validation summaries, and reproducibility metadata may be committed.

A separate long-lived `develop` branch is intentionally not used.

## 2. Branch naming

Use one branch for one scientific or engineering idea.

- `fix/<topic>` — bug or correctness fix
- `perf/<topic>` — performance-only optimization
- `physics/<topic>` — physical model or definition change
- `analysis/<topic>` — derived observable / analysis change
- `validation/<topic>` — validation protocol or convergence tooling
- `docs/<topic>` — documentation only
- `chore/<topic>` — repository maintenance

Do not mix an unrelated cleanup, physics change, and performance optimization in one PR.

## 3. Change classes and required evidence

### Documentation only

Examples: prose, comments, workflow documentation.

Required:
- no source or numerical behavior change;
- review links and commands for correctness when applicable.

### Engineering / performance

Examples: logging, caching, data layout, faster but mathematically equivalent evaluation.

Required:
- relevant unit/regression tests;
- full Release CTest before merge;
- before/after benchmark when performance is claimed;
- a check that the optimized path remains statistically / numerically consistent with the reference path when floating-point execution can diverge.

Performance diagnostics are not physics uncertainties unless the estimator changes.

### Numerical method

Examples: interpolation, tolerances, quadrature, binning, integrator behavior.

Required:
- full Release CTest;
- targeted deterministic regression tests;
- a convergence or reference comparison appropriate to the changed numerical method;
- numerical settings written to metadata when they can affect results.

Do not promote a faster setting to production only because two noisy Monte Carlo runs are statistically compatible.

### Physics definition

Examples: capture definition, orbit removal, sampling law, scattering model, annihilation treatment, weighted histories.

Required:
- explicit statement of the old and new physical definitions;
- updated documentation and output/schema provenance;
- deterministic tests for invariants where possible;
- new physics-facing validation;
- explicit statement that old production results are or are not compatible.

A physics-definition change normally invalidates previous production results unless compatibility is demonstrated.

### Analysis / observable

Examples: annihilation integrals, gamma-ray projection, neutrino likelihoods, uncertainty estimators.

Required:
- fixture tests and mathematical closure checks;
- backward compatibility decision for stored products;
- independent validation of nonlinear estimators when finite-sample bias can matter.

## 4. Pull-request workflow

1. Synchronize from `main`.
2. Create a short-lived branch.
3. Make the smallest coherent change.
4. Run focused tests during development.
5. Run the full required gate before declaring the PR ready.
6. Push the branch and open a PR using the repository template.
7. Let CI pass.
8. Prefer squash merge so `main` contains one coherent commit per PR.
9. Delete the short-lived branch after merge.

Direct pushes, force pushes, and history rewriting on `main` are not part of the normal workflow.

For a single-maintainer repository, external approval need not be mandatory, but the PR and CI gate should still be used as the permanent scientific change record.

## 5. Automated CI versus scientific validation

CI answers: **did the code break?**

The CI gate runs a Release build and the tracked CTest suite, including tests that are available on the runner. It is expected to cover unit tests, analysis contracts, and supported MPI regressions.

Scientific validation answers: **is the result trustworthy at the required precision?**

Examples include:

- sample-size convergence;
- independent-seed studies;
- numerical-tolerance convergence;
- direct versus interpolated scattering-rate comparisons;
- radial-resolution / rebin studies;
- outer-boundary sensitivity;
- thermal-limit recovery.

Scientific validation is not required on every small commit, but it is required before freezing a production version affected by the relevant method.

Keep protocols and compact summaries under `validation/`. Keep raw run directories outside Git.

## 6. Current transport invariants

These are current model definitions, not universal truths. Changing any of them is a `physics/` or `numerical/` task and must be explicit.

- incident sampling reference: 1100 AU;
- deterministic injection checkpoint: 2 solar radii;
- numerical/analytic matching surface: 1 solar radius;
- default bound-orbit removal: 1100 solar radii;
- native radial grid: 0.001 solar-radius bins through 1.1 solar radii, then 2% geometric width growth capped at 10 solar radii;
- capture normalization comes from an independent fixed-incident run;
- complete captured transport uses physical escape or explicit outer-orbit removal as complete outcomes;
- computational prefixes do not normalize a production occupation;
- outer-orbit removal is not evaporation;
- no finite solar-age cutoff is part of the current complete-history definition.

The README and implementation are authoritative if this list becomes stale; update this file in the same PR as any intentional definition change.

## 7. Production rules

A paper-quality production campaign should satisfy all of the following:

- exact commit or release tag recorded;
- clean source tree;
- Release build and build configuration recorded;
- source hash and numerical settings recorded in output metadata;
- explicit nonzero seeds and MPI rank count recorded;
- production acceptance gates satisfied;
- no unresolved captured-transport numerical or computational failures;
- required scientific convergence studies completed for the chosen setup.

Do not use `git pull && run` as a production procedure. Check out an exact tag or commit, build it, run the test gate, then launch the campaign.

Recommended tag convention:

- `vX.Y.Z` for software releases;
- an additional descriptive tag such as `paper-v1` may point to the exact release/commit used for a manuscript production campaign.

## 8. Data and validation artifacts

Do not commit raw Monte Carlo run directories, large generated tables, scheduler logs, or local binaries.

For an important production/validation campaign, preserve a compact manifest containing at least:

- campaign name;
- code commit / tag;
- source hash;
- physical parameter grid;
- numerical settings;
- seeds;
- MPI layout;
- analysis version;
- location / identifier of external raw data;
- completion and validation status.

A result is not considered reproducible from a plot alone.

## 9. GitHub repository settings

For `main`, enable a branch ruleset / protection rule with:

- pull request required before merge;
- the `build-test` CI job required;
- force pushes blocked;
- branch deletion blocked.

For a single maintainer, a mandatory external approval is optional.

Use squash merge as the default merge method for development PRs.

## 10. Before a production freeze

Before creating a production tag:

- [ ] full Release CTest passes;
- [ ] relevant scientific validation is current for the exact numerical setup;
- [ ] known limitations are written down;
- [ ] configuration and metadata fields fully identify the run;
- [ ] the source tree is clean;
- [ ] the release / campaign manifest is archived;
- [ ] the tag points to the exact reviewed commit.
