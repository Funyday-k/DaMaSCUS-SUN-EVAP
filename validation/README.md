# Physics validation workflow

This directory complements the ordinary unit tests with physics-facing
verification and production-scale convergence checks.

## Fast gate

The `test_Physics_Validation` CTest target checks:

- fourth-order convergence of the RKF45 propagator against an analytic
  eccentric Kepler orbit;
- monotonic/positive structure of the imported solar profiles;
- the area law for sampled incident impact parameters;
- statistical agreement of dark-photon angle samples with their CDF;
- the analytic Maxwell mean-speed limit for a stationary DM particle.

Run it with:

```bash
cmake --build build --target test_Physics_Validation --parallel
ctest --test-dir build --output-on-failure -L physics-validation
```

## Production convergence matrix

`physics_validation.py` runs the same physical configuration for a matrix of
scattering-rate interpolation grids and PRNG seeds. Grid `0` is the direct-rate
baseline. Each interpolated run is checked against the same-seed baseline.

```bash
python3 validation/physics_validation.py \
  --executable build/src/DaMaSCUS-SUN \
  --config /path/to/physics.cfg \
  --output validation-results/m0p5-s36 \
  --grids 0 1000 2000 \
  --seeds 271828 314159 161803 \
  --sample-size 1000 \
  --ranks 8
```

The source configuration is never modified. For every matrix point the tool
writes an isolated config, log, and result directory. The final
`physics_validation_report.json` records all metrics, checks, thresholds, and
the overall pass/fail decision.

Default gates are:

- no early stop and numerical-failure fraction at most `1e-3`;
- overlapping 95% Wilson intervals for the valid capture rate;
- average scattering counts within 10%;
- captured radial distributions within total-variation distance 0.10;
- complete-evaporation fractions compatible within a two-proportion
  three-sigma test;
- observed evaporation-time medians within 20% and 90th percentiles within
  30% when both runs contain at least 30 complete events.

Use `--require-evaporation-events` when evaporation is the primary observable;
otherwise insufficient event counts are explicitly reported as
`not_evaluated`, not silently treated as evidence. Use `--dry-run` to inspect
the generated matrix configurations before launching expensive jobs, and
`--reuse-existing` to regenerate a report from completed matrix directories.

These thresholds are validation defaults, not universal physics constants.
Tighten them when the Monte Carlo sample is large enough, and record any
project-specific threshold changes with the resulting report.
