# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun, with capture- and evaporation-focused
extensions.

The format-2 `bincount.tsv` transport contract, configuration and validation scope are documented below.

## Overview

DaMaSCUS-SUN-EVAP builds on
[DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN) and keeps its central
Monte Carlo picture: dark matter particles are propagated through the solar
potential, scatter on solar targets, and are classified from their trajectory
history. This branch is organized around production workflows for low-mass dark
matter capture and evaporation studies rather than broad direct-detection scans.

The current code path is centered on two practical modes:

- **Capture mode**: a fast capture-rate workflow. It terminates a trajectory once
  a post-scatter bound state is identified and avoids the full evaporation and
  histogram output path.
- **Parameter-point simulation**: the main evaporation workflow for one mass and
  cross section. It accumulates time-weighted radial histograms, records complete
  captured histories ending in escape or outer-orbit removal, and can emit wall-clock snapshot progress files for
  long MPI runs.

The older parameter-scan machinery is still present, but the most actively
maintained outputs in this branch are the capture summary and the single
parameter-point evaporation products.

## Main Changes

Compared with the upstream DaMaSCUS-SUN workflow, this branch emphasizes:

- in-memory radial bincount accumulation instead of writing full trajectory
  files;
- capture detection from the first post-scatter negative-energy state;
- one final radial sufficient-statistics file for scientific transport;
- optional richer survival diagnostics behind explicit diagnostic paths;
- MPI-aware snapshot output for long parameter-point jobs;
- safeguards for pathological trajectories and low-capture-rate runs;
- server-friendly local configuration conventions, with generated binaries,
  job scripts, and run configs kept outside version control under `bin/`.

## Build And Deployment

### Dependencies

- CMake 3.12 or newer and Git
- C++14-capable compiler
- OpenMPI or MPICH, including development headers
- Boost 1.65 or newer
- libconfig++ development headers and library
- Python 3 only when building/running the full test suite

On Ubuntu/Debian, install the native prerequisites with:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git libboost-dev libconfig++-dev libopenmpi-dev openmpi-bin python3
```

On macOS with Homebrew:

```bash
brew install cmake boost libconfig open-mpi
```

The first CMake configure downloads `obscura` v1.0.1 and `libphysica` v0.1.2
under the build tree's `_deps/` directory. Their exact Git commits are pinned by the root CMake project.
GoogleTest is fetched only when `BUILD_TESTING=ON`. Internet access is therefore
needed on the first configure unless those FetchContent source directories have
already been populated.

### Local Build

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/install" \
  -DBUILD_TESTING=OFF \
  -DCODE_COVERAGE=OFF
cmake --build build --config Release --parallel
cmake --install build --config Release
```

This produces a self-contained project install layout (apart from the native
MPI and libconfig++ shared libraries):

```text
install/
├── bin/DaMaSCUS-SUN
└── share/DaMaSCUS-SUN/
    ├── model_agss09.dat
    └── examples/quickstart.cfg
```

Verify the fresh install with the tracked, intentionally tiny configuration:

```bash
mpirun -np 1 ./install/bin/DaMaSCUS-SUN \
  ./install/share/DaMaSCUS-SUN/examples/quickstart.cfg
```

The installed program locates `model_agss09.dat` relative to its own executable,
not relative to the checkout or the current working directory. The complete
`install/` directory can therefore be moved to another location on the same
machine. Do not copy only the executable; copy the whole install prefix.

Native binaries still use the MPI and libconfig++ libraries from the build
machine. When moving to a machine with a different OS, CPU architecture, MPI
implementation, or module stack, rebuild there. On module-based clusters, load
the same compiler/MPI/libconfig modules for both build and execution.

### Cluster Deployment

A typical cluster checkout follows the same pattern:

```bash
git clone https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP.git
cd DaMaSCUS-SUN-EVAP
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/install" \
  -DBUILD_TESTING=OFF \
  -DCODE_COVERAGE=OFF
cmake --build build --config Release --parallel
cmake --install build --config Release
```

Keep machine-specific configuration, batch scripts, and scheduler output outside
the install prefix (the repository's ignored `bin/` directory is one option).

Run a parameter point with MPI:

```bash
mpirun -np 8 ./install/bin/DaMaSCUS-SUN /absolute/path/to/config.cfg
```

For a scheduler, wrap the same executable/config pair in the local batch script
used by that machine.

If a nonstandard packaging layout stores the solar table elsewhere, set either
an exact file or its containing directory:

```bash
export DAMASCUS_SUN_SOLAR_MODEL=/absolute/path/model_agss09.dat
# or: export DAMASCUS_SUN_DATA_DIR=/absolute/path/to/data
```

To build and run the test suite, configure separately with
`-DBUILD_TESTING=ON`, build, then run `ctest --test-dir build --output-on-failure`.

Incident particles are sampled on a fixed 1100 AU reference sphere, mapped
analytically through 2 R_sun to the 1 R_sun solar matching surface, and then
propagated through the Sun. Bound exterior orbits are removed on their first
outward crossing of the independently configured 1100 R_sun cutoff. Incoming
shell moments are recorded from the diagnostic reference (up to 1100 R_sun)
to the solar surface; capture-conditioned residence starts only after capture.
Snapshot files are progress products.

## Configuration

Configuration files use libconfig syntax. The most important controls are:

| Setting | Meaning |
| --- | --- |
| `run_mode` | `"Parameter point"` for the main evaporation workflow, `"Capture"` for capture-rate runs, or `"Parameter scan"` for detector-limit scans. |
| `sample_size` | In Parameter point mode, the exact number of complete captured histories (escape or outer removal). In Capture mode, the exact number of incident trials. Failed histories are recorded and skipped; Transport continues until the complete-capture target or the explicit attempt budget is reached. |
| `outer_removal_radius_rsun` | Bound-orbit removal radius in R_sun, default 1100; must exceed the native 1.1 R_sun grid. The obsolete `outer_boundary_radius_au` key is rejected. |
| `fixed_seed` | Optional non-negative PRNG seed. `0` or an omitted setting uses nondeterministic seeding; a nonzero value is expanded independently by MPI rank. |
| `max_trajectories` | Optional hard cap on generated trajectories. `0` or unset means no trajectory-count cap. |
| `interpolation_points` | Legacy square scattering-rate grid size. It remains supported; the three `rate_*` settings below override its corresponding defaults. `0` disables interpolation when no rectangular-grid settings are supplied. |
| `rate_radius_points` | Optional number of radial rate-grid points. Defaults to `interpolation_points`. Explicit rate grids must use `(0, 0)` to disable interpolation or at least two points in both dimensions. |
| `rate_speed_points` | Optional number of speed rate-grid points. Defaults to `interpolation_points`; it follows the same explicit-grid validation as `rate_radius_points`. |
| `rate_max_speed` | Optional maximum tabulated DM speed in natural units (`0.02` means `0.02c`), constrained to `(0, 0.75]` and defaulting to `0.75`. Faster queries fall back to the direct rate and are counted. When the grid is disabled, this input is ignored and the bincount header records the actual table limit `0`. |
| `output_dir` | Root directory for generated result folders; a trailing `/` is optional. A relative path is resolved from the process working directory, so production batch jobs should normally use an absolute path. |
| `DM_mass` | Dark matter mass in GeV. |
| `DM_cross_section_nucleon` | DM-nucleon cross section in cm^2. |
| `DM_cross_section_electron` | DM-electron cross section in cm^2 where relevant. |
| `maximum_number_of_scatterings` | Per-trajectory computational cutoff. Cutoff-terminated captures are not treated as clean physical evaporation events. |
| `snapshot_enabled` | Enables intermediate wall-clock progress reports for parameter-point runs. Disabled automatically in capture mode. |
| `snapshot_interval` | Positive integer wall-clock spacing, in seconds, for snapshot reports. Defaults to 60 seconds when snapshots are enabled. |
| `max_trajectory_wall_time_sec` | Optional per-trajectory wall-time guard. Snapshot recorder overhead is excluded from this budget. |

MPI trajectory scheduling uses a dynamic RMA work queue. Every rank claims one
trajectory at a time and releases its slot immediately on completion, so fast
ranks continue working without waiting at a per-batch collective for the
slowest trajectory. The queue maintains
`accepted_samples + in_flight <= sample_size`; therefore the final exact-target
tail may temporarily leave excess ranks idle, but
`capture_target_overshoot` remains zero. The terminal summary reports work claims and peak in-flight counts.

Rank 0 also advances MPI during trajectory propagation with a nonblocking
`MPI_Iprobe`, at most once per millisecond. This is needed on MPI transports
that do not progress passive-target RMA while the window owner is computing:
otherwise rank 0's first long trajectory can stall every other rank's first
claim. Progress runs on the main thread, independently of snapshots, and all
ranks still compute trajectories. The scheduling algorithm is unchanged by the output contract.

The exact-target rule still applies: `sample_size = 1` allows only one active
trajectory, and fewer than 32 remaining target slots cannot keep 32 ranks busy.
Increasing MPI ranks alone does not remove that intentional tail limit.

For reproducible MPI runs, a nonzero fixed seed is expanded by rank as
`base_seed + 1000003 * mpi_rank`. Ordinary Capture and Parameter point runs
record numerical failures and computational truncations and continue issuing
work. Initial-shift failures also continue; no failure-fraction threshold stops
the simulation. A failed Transport history contributes neither residence nor
complete-path moments. Transport finishes after `sample_size` complete captured
histories; Capture still finishes after exactly `sample_size` incident attempts,
including failed attempts. `max_trajectories` remains an explicit total-attempt
budget. Without that budget, Transport continues seeking complete captures.

`target_reached` reports completion only. There is no zero-error acceptance gate.
Exit code 2 means the configured attempt budget prevented reaching the target;
filesystem/configuration errors still fail normally. Rank assignment can vary
with MPI scheduling; bitwise comparisons should use a single rank and no wall-time cutoff.

Remove `production_mode`, `thermal_validation_mode`, `trajectory_summary_enabled`,
`trajectory_events_enabled`, `trajectory_trace_rate` and `trajectory_trace_seed`
from ordinary configuration files. These obsolete keys now produce an explicit
migration error. A pilot uses the same scientific contract with a smaller
`sample_size`; it does not enable diagnostics or change error handling.

A trajectory can become physically bound only at a scattering. If an
uncaptured trajectory acquires negative energy during scatter-free propagation,
the run classifies that trajectory as a numerical failure instead of allowing
repeated bound Kepler returns to stall its MPI batch.

## Outputs

Capture (`run_mode = "Capture"`) is a fixed **incident-count** normalization run.
Stdout contains exactly one `CAPTURE_RESULT_JSON={...}` line (capture schema 2);
human-readable logs go to stderr. No result directory, snapshot, copied cfg or
diagnostic file is created, including on a failed run. The record retains
target completion, physics/numerics, seed, MPI ranks, failure counts,
`N_inj`, `N_capt`, `f_cap`, `C_geom_s_inv`, `C_capture_s_inv` and 64 count blocks.
Failed histories do not stop the fixed-injection run. `N_valid` and
`N_unclassified` identify classified and failed attempts. For compatibility with
the fixed-incident estimator, `f_cap` remains `N_capt / N_inj`, explicitly labelled
`f_cap_denominator = all_injected_trials`; failure counts are retained. The
`target_reached` field replaces the old `production_accepted` veto.

Ordinary Transport (`run_mode = "Parameter point"`) has one output contract:

```text
results_<log10_mass_GeV>_<log10_sigma_p_cm2>/
├── bincount.tsv
└── snapshot/
```

The result directory must be new or empty. A nonempty result directory is
rejected before simulation, without changing its files. Use a new `output_dir`
for each run, including retries and concurrent jobs. This protects old snapshots
and prevents an earlier accepted file from being mistaken for a new result.
There is no automatic cleanup or migration of historical result directories.

`bincount.tsv` is written to a temporary file, checked, closed, and atomically
renamed. Its header contains `# key = value` entries for effective physics and
numerics, actual seed, MPI ranks, global counts, target completion and stopping reason.
Cross sections explicitly distinguish proton, neutron and electron values.
The input model settings retain their libconfig units, with halo velocities
labelled km/s. No source fingerprint or Git compatibility gate is required.

The same file contains exactly 64 comment records:

```text
# block_count_columns = block N_injected N_ever_captured N_never_captured N_residence_samples
# block_count = 0 3041 157 2883 156
```

The numeric table has eleven columns, named in its `# columns = ...` comment:

| Columns | Meaning and units |
| --- | --- |
| `block`, `bin` | Zero-based block and radial-bin indices. |
| `r_low_rsun`, `r_high_rsun` | Shell boundaries in R_sun. |
| `captured_residence_dt_sum_s` | Sum of per-history residence after first capture [s]. |
| `captured_residence_dt_sq_sum_s2` | Sum of squared per-history captured residence [s²]. |
| `captured_residence_v2dt_sum_km2_s` | Captured-residence integral of speed squared [km²/s]. |
| `ever_captured_path_dt_sum_s`, `ever_captured_path_dt_sq_sum_s2` | First and second time sums for complete recorded ever-captured paths [s, s²]. |
| `never_captured_path_dt_sum_s`, `never_captured_path_dt_sq_sum_s2` | First and second time sums for complete recorded never-captured paths [s, s²]. |

Rows are ordered by block, then bin, including zero rows. The existing grid is
unchanged: 0.001 R_sun shells through 1.1 R_sun, then shell widths grow by 2%,
capped at 10 R_sun and clipped at the removal boundary. The default 1100 R_sun
cutoff gives 1626 bins and 104064 rows. For example:

```python
import numpy as np
radial = np.loadtxt("bincount.tsv", comments="#")
# columns: block, bin, r_low, r_high, res_dt, res_dt2, res_v2dt,
#          ever_dt, ever_dt2, never_dt, never_dt2
```

Captured residence ends at validated escape on the solar matching surface or
outer removal. Ever-captured paths include inbound, pre-capture, residence and
outgoing components; these are added **before squaring each history's bin time**.
Never-captured paths include scattered and unscattered completed escapes.
Incoming/outgoing paths are recorded within `R_path_reference_rsun`; captured
bound paths extend to `R_remove_rsun`. This preserves the existing boundaries.

Counts distinguish all ever-captured histories from complete captured histories:
`N_ever_captured` includes captures that later failed, while `N_residence_samples`
is the denominator for both captured-residence and ever-captured complete-path
moments. Each block includes that denominator as its fourth count.
`N_never_captured` counts completed never-captured histories only.

For every global/block population,
`N_residence_samples + N_never_captured + excluded = N_injected`.
The header records `N_excluded_trajectories`; a block's excluded count is the
injected count minus its two complete-history counts. `N_unclassified` records
failed histories that were never classified as captured. Captured histories that
later failed are the difference `N_ever_captured - N_residence_samples`.
When there are no failures, these reduce to the original two-population closure.
The eleven radial columns are unchanged. Format 2 distinguishes this count
contract and `target_reached` from the former zero-error `accepted` field.

A run that reaches its target writes `target_reached = true` even when failure
counts are nonzero. Exhausting an explicit attempt budget writes
`target_reached = false` and exits 2. An I/O failure exits nonzero collectively
and does not publish a partial file.

Local analysis uses the **independent fixed-injection Capture** probability and
conditional Transport moments. Do not estimate population fractions from the
fixed-captured-count Transport ratio. The time-square sums support per-bin pair
estimators; for a population with N >= 2, use `(S1*S1 - S2)/(N*(N-1))`.
Block deletion retains correlated radial first moments for jackknife analysis.
Per-history cross-bin second moments and full velocity distributions are not
part of this contract.

Analysis and plotting live in the sibling DaMaSCUS-SUN repository. Readers of
old `metadata.json`/`radial_blocks.tsv` products must migrate to format 2;
there is no dual-write compatibility mode. Keep cfg files and scheduler logs
outside the result directory.

Local diagnostics and thermal shape validation remain explicit tools:

```bash
DaMaSCUS-SUN CONFIG.cfg --diagnostic
DaMaSCUS-SUN CONFIG.cfg --thermal-validation
replay-invalid-trajectory CONFIG.cfg invalid_trajectories.tsv RANK TRAJECTORY_ID
```

Both local entry points require `run_mode = "Parameter point"` and write under
`output_dir/diagnostics/results_.../`. They preserve legacy diagnostic/replay
reports, including `bincount.txt`, `evaporation_times.txt`,
`residence_jackknife_blocks.tsv`, `run_metadata.json`,
`diagnostic_trajectory_summary.tsv`, `trajectory_events.tsv` and
`invalid_trajectories.tsv`. These are separate local workflows; they do not
publish a scientific `bincount.tsv`. The diagnostic CLI traces all histories;
tests can configure narrower tracing through the internal API. Thermal shape
validation may retain computationally truncated residence prefixes and is
never accepted for absolute production analysis.

The replay helper restores the recorded RNG state and shifted initial condition
and prints the original/replayed termination details. No trajectory ledger or
replay data is collected by ordinary scientific runs.

When snapshots are enabled, intermediate files are written under `snapshot/`:

- `snapshot_{time}s.txt`: cumulative progress report at the snapshot wall time.
  Its commented `[MPI rank status]` table reports each rank's activity, local
  trajectory ID, trajectory wall time, simulated elapsed time, scattering
  count, and the rank-local observation time. An in-progress trajectory that
  has already captured contributes its accumulated residence prefix
  provisionally. After each bound exterior arc returns, a forced publication
  includes the complete round-trip residence integral, extending the histogram
  to the orbit's apoapsis. Status rows remain comments so data readers see only
  bincount bins.
- `snapshot_{time}s_evaporation_times.txt`: complete valid evaporation events
  first published by that checkpoint, sorted by `lifetime_unbinding_sec`.
  An event committed concurrently with a snapshot boundary is assigned once to
  the next checkpoint rather than being dropped.

Snapshot files are progress diagnostics. They do not replace the final
post-reduction `bincount.tsv`, and they are
not restart checkpoints. A report can temporarily have `snapshot_status =
partial` while ranks publish their state. If a rank misses the deadline, the
report remains incomplete rather than reconstructing state that was not
captured. Final rank states are retained when the last merge is incomplete.
Snapshot reports likewise expose attempted, classified, unresolved, and
`numerical_failures` counts plus raw and valid capture-rate intervals. The
uncaptured-bound free-flight guard records its failures in these cumulative
reports without printing one warning for every failed trajectory.

The executable requests `MPI_THREAD_FUNNELED`; only the main thread calls MPI,
while one heartbeat thread per rank performs local state copies and file I/O.
Concurrent jobs must use distinct output directories.

Snapshot checkpoint I/O is supported on homogeneous POSIX (Linux/macOS) MPI
nodes sharing the output filesystem. The binary checkpoint representation is
local to a run and is not a portable interchange or restart format.

The `snapshot/` directory remains after completed or budget-limited Transport runs,
even when empty because snapshots were disabled or the run ended before its
first interval. Existing snapshot contents, heartbeat behavior and checkpoint
formats are unchanged. Capture never creates this directory.

## Numerical tests

The ordinary test suite now includes a `physics-validation` gate for analytic
Kepler convergence, solar-profile structure, incident impact-parameter
sampling, scattering-angle sampling, and the Maxwell thermal-speed limit:

```bash
ctest --test-dir build --output-on-failure -L physics-validation
```

## Citation

If this branch is used in analysis, cite this repository using `CITATION.cff`
and cite the original DaMaSCUS-SUN work where appropriate.

Useful references:

- T. Emken and C. Kouvaris,
  [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN)
- Garani and Palomares-Ruiz, evaporation of dark matter in the Sun

## License

This project is distributed under the MIT License. See `LICENSE` for the
upstream and modification copyright notices.
