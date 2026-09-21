# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun, with capture- and evaporation-focused
extensions.

The current schema-10 transport definitions, configuration, output contract, and validation scope are documented below.

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
- compact final evaporation-time output for complete valid unbinding events;
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
| `sample_size` | In Parameter point mode, the exact number of complete captured histories (escape or outer removal). In Capture mode, the exact number of incident trials. Invalid histories invalidate production even if replacements reach the target. |
| `production_mode` | Stop issuing new trajectories after any failed/truncated history, drain in-flight work, and exit nonzero; an unmet target also fails. Capture stdout JSON or transport metadata records acceptance; diagnostics/replay are disabled. |
| `outer_removal_radius_rsun` | Bound-orbit removal radius in R_sun, default 1100; must exceed the native 1.1 R_sun grid. The obsolete `outer_boundary_radius_au` key is rejected. |
| `thermal_validation_mode` | Separate Parameter point shape workflow allowing computationally limited histories; never absolute production. |
| `fixed_seed` | Optional non-negative PRNG seed. `0` or an omitted setting uses nondeterministic seeding; a nonzero value is expanded independently by MPI rank. |
| `max_trajectories` | Optional hard cap on generated trajectories. `0` or unset means no trajectory-count cap. |
| `interpolation_points` | Legacy square scattering-rate grid size. It remains supported; the three `rate_*` settings below override its corresponding defaults. `0` disables interpolation when no rectangular-grid settings are supplied. |
| `rate_radius_points` | Optional number of radial rate-grid points. Defaults to `interpolation_points`. Explicit rate grids must use `(0, 0)` to disable interpolation or at least two points in both dimensions. |
| `rate_speed_points` | Optional number of speed rate-grid points. Defaults to `interpolation_points`; it follows the same explicit-grid validation as `rate_radius_points`. |
| `rate_max_speed` | Optional maximum tabulated DM speed in natural units (`0.02` means `0.02c`), constrained to `(0, 0.75]` and defaulting to `0.75`. Faster queries fall back to the direct rate and are counted. When the grid is disabled, this input is ignored and metadata records the actual table limit `0`. |
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
`capture_target_overshoot` remains zero. Output headers report
`mpi_scheduler_work_claims` and `mpi_scheduler_peak_in_flight`.

Rank 0 also advances MPI during trajectory propagation with a nonblocking
`MPI_Iprobe`, at most once per millisecond. This is needed on MPI transports
that do not progress passive-target RMA while the window owner is computing:
otherwise rank 0's first long trajectory can stall every other rank's first
claim. Progress runs on the main thread, independently of snapshots, and all
ranks still compute trajectories. Updated final headers identify this path as
`mpi_scheduler_progress = main_thread_iprobe_v1`.

The exact-target rule still applies: `sample_size = 1` allows only one active
trajectory, and fewer than 32 remaining target slots cannot keep 32 ranks busy.
Increasing MPI ranks alone does not remove that intentional tail limit.

For reproducible MPI runs, a nonzero fixed seed is expanded by rank as
`base_seed + 1000003 * mpi_rank`. Computational cutoffs are tracked separately
from physical right-censoring so that final evaporation-time files contain only
complete valid unbinding events. Numerical failures and computational
truncations are recorded and replaced; their accumulated fraction does not
stop the work queue. Use `max_trajectories` when an explicit attempt budget is
required.

A trajectory can become physically bound only at a scattering. If an
uncaptured trajectory acquires negative energy during scatter-free propagation,
the run classifies that trajectory as a numerical failure instead of allowing
repeated bound Kepler returns to stall its MPI batch.

## Outputs

Capture (`run_mode = "Capture"`) is a fixed **incident-count** normalization run.
It never creates a result directory, snapshot, copied cfg or diagnostic file.
Rank zero prints exactly one `CAPTURE_RESULT_JSON={...}` line (capture schema 1),
including accepted/rejected status, effective physics/numerics, seed and MPI ranks,
failure counts, `N_inj`, `N_capt`, `f_cap`, `C_geom_s_inv`, `C_capture_s_inv`, and
64 `[N_inj, N_capt]` blocks. Failed/incomplete runs still print the record and exit 2.
Keep the external cfg and Slurm stdout.

Production transport (`run_mode = "Parameter point"`, `production_mode = true`)
writes only `metadata.json` (schema 10), `radial_blocks.tsv`, `block_counts.tsv`,
`trajectory_summary.tsv`, `orbit_class_blocks.tsv`, `incident_inbound.tsv`,
`termination_counts.tsv`, `solar_reference.tsv`, and the enabled `snapshot/`.
`trajectory_summary.tsv` is the reduced scientific summary, independent of the
optional diagnostic summary. Metadata is the canonical effective configuration,
including maximum scatterings, trajectory wall-time budget, production and thermal
flags. No `input.cfg` is copied. Snapshots report progress; `restart_supported=false`.
Runs shorter than `snapshot_interval` may leave an empty snapshot directory.

Analyze locally with an independent capture log or an extracted JSON record:

```bash
python scripts/analyze_point.py TRANSPORT_DIR --capture-log logs/capture.out
# Alternatively: --capture-json capture_result.json
python scripts/prepare_transport_runs.py campaign --phase pilot
```

The generator writes `configs/POINT/seedN/{capture,transport}.cfg` and a manifest
with log paths under `logs/` and transport destinations under `results/`.
Capture disables snapshots; transport enables them. Both production configurations
disable trajectory diagnostics. No jobs are launched by the generator.
The analyzer requires schema 10 plus capture schema 1, checks physics, solar-model
identifier, numerical settings and disjoint rank seeds, and combines independent
capture/transport delete-block jackknife uncertainties. Hashes are not validated.
Analyze older schema 8/9 data using the code revision that produced it.

Local transport (`production_mode = false`) also supports the legacy products below.
Enable `trajectory_summary_enabled`, `trajectory_events_enabled`, and
`trajectory_trace_rate` as needed for diagnosis and replay:

- `bincount.txt`: legacy capture-conditioned residence and velocity-moment output.
  The grid is uniform at 0.001 R_sun through 1.1 R_sun, then grows by 2% per shell
  with a 10 R_sun width cap, clipped at the removal surface. Analytic exterior
  arcs use a round trip or a one-way removal arc as appropriate. Computational
  and numerical failures do not enter production residence. Prefixes may appear
  only in the explicitly labelled thermal shape workflow. Use schema-10 products
  for the new independent capture normalization.
- `evaporation_times.txt`: compact complete-event table with
  `rank trajectory_id lifetime_unbinding_sec r_capture_Rsun E_capture_eV
  dE_capture_eV`, followed by the number of negative-energy exterior arcs,
  the first/last/maximum osculating Kepler periods at outward `1 R_sun`
  crossings, and the corresponding first/last/maximum analytic exterior
  return times. It is sorted by
  `lifetime_unbinding_sec` with `rank trajectory_id` tie-breakers.
- `residence_jackknife_blocks.tsv`: exactly 64 deterministic blocks assigned by
  `splitmix64(base_seed, rank, trajectory_id) % 64`. Each block contains
  attempted, captured, completed uncaptured escape, accepted residence,
  invalid, and outer-orbit-removal counts
  plus its full radial `dt` and `v^2 dt` histograms. The writer refuses to publish the file unless every
  scalar count and every radial bin closes against `bincount.txt`. This legacy joint-run product does not replace the independent capture and
  transport blocks used by the schema-10 analysis.
- `invalid_trajectories.tsv`: local-only replayable ledger for trajectories
  excluded by numerical or computational validity rules. It is header-only
  when no invalid trajectory occurred. Each row records the failure stage,
  exact termination reason and numerical-failure detail, boundary/reference
  energy diagnostics, capture/survival state, final kinematics, shifted initial
  condition, and the `std::mt19937` states before initial-condition generation
  and before trajectory simulation.

Replay one ledger row with the installed helper:

```bash
replay-invalid-trajectory CONFIG.cfg invalid_trajectories.tsv RANK TRAJECTORY_ID
```

The helper restores the recorded pre-simulation RNG state and shifted initial
condition, reruns the current trajectory implementation, and prints the
original/replayed reason, failure detail, final state, and diagnostic-event
count. It reads only the current ledger schema; old ledger compatibility is not
provided.

The terminal summary and `bincount.txt` header report captured and uncaptured
counts for every `TrajectoryTerminationReason` and every concrete numerical
failure detail. This breakdown is always available; it is not gated by
`trajectory_summary_enabled`.

Local trajectory diagnostics are enabled with `trajectory_summary_enabled = true`.
This adds `run_metadata.json`, `diagnostic_trajectory_summary.tsv`, and
`trajectory_events.tsv` without changing the capture or evaporation state
definitions. Set `trajectory_events_enabled = true` and
`trajectory_trace_rate` in `[0, 1]` to select lifecycle traces by a stable hash;
set `trajectory_trace_seed` to keep that selection identical across runs. The
selection does not consume the physics RNG. Traced trajectories include the
complete pre-initial-condition `std::mt19937` state and shifted initial
condition, plus real-time scatter, state-transition, solar-crossing, escape,
censoring, and numerical-failure events.

The bound-exit period is the point-mass osculating Kepler period inferred from
the negative-energy state at the outward `1 R_sun` matching surface. The
exterior elapsed time is the physically used analytic travel time: through
apoapsis to the inbound matching surface for every bound exterior arc. These
are kept separate because the osculating full period includes a point-mass continuation through
the solar interior, whereas the simulation uses the extended solar potential
there.

The `bincount.txt` and snapshot report headers expose both `capture_rate_raw`
(captured / all attempted) and `capture_rate_valid` (captured / physically
classified), with separate standard errors and Wilson intervals.

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
post-reduction `bincount.txt` and `evaporation_times.txt` products, and they are
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

Capture-mode runs skip the full output path and print the capture summary
instead.

## Physics Validation

The ordinary test suite now includes a `physics-validation` gate for analytic
Kepler convergence, solar-profile structure, incident impact-parameter
sampling, scattering-angle sampling, and the Maxwell thermal-speed limit:

```bash
ctest --test-dir build --output-on-failure -L physics-validation
```

Production results should additionally pass the slower grid-and-seed
convergence matrix described in
[`validation/README.md`](validation/README.md). It compares direct scattering
rates with interpolated grids using capture-rate confidence intervals, average
scattering counts, radial distributions, evaporation-time statistics, and the
numerical-failure fraction.

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
