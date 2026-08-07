# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun, with capture- and evaporation-focused
extensions.

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
  valid evaporation events, and can emit wall-clock snapshot progress files for
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

The radial residence domain is a compile-time setting because its bin count is
part of the fixed histogram layout. Keep alternate cutoffs in separate build
directories and provide both the cutoff and the matching number of exterior
shells. The default is `1 AU` with `423` exterior bins. For example, a
5.2-AU comparison build uses:

```bash
cmake -S . -B build-5p2au \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DDAMASCUS_RADIAL_DOMAIN_MAX_AU=5.2 \
  -DDAMASCUS_EXTERIOR_BINS=527
cmake --build build-5p2au --target DaMaSCUS-SUN --parallel
```

Grid initialization rejects a mismatched bin count when the generated shells
either stop short of the chosen domain or reach it before the final bin.

## Configuration

Configuration files use libconfig syntax. The most important controls are:

| Setting | Meaning |
| --- | --- |
| `run_mode` | `"Parameter point"` for the main evaporation workflow, `"Capture"` for capture-rate runs, or `"Parameter scan"` for detector-limit scans. |
| `sample_size` | In normal mode, the exact target number of complete, valid evaporation events whose bound exterior orbit stays within the default 1-AU radial domain. Captures physically removed at 1 AU contribute residence statistics but not the evaporation-event target. Invalid captures are replaced. In capture mode, this is the exact target number of captures. |
| `fixed_seed` | Optional non-negative PRNG seed. `0` or an omitted setting uses nondeterministic seeding; a nonzero value is expanded independently by MPI rank. |
| `max_trajectories` | Optional hard cap on generated trajectories. `0` or unset means no trajectory-count cap. |
| `interpolation_points` | Scattering-rate interpolation grid size. `0` disables interpolation; production runs should compare representative values before fixing this. |
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

For reproducible MPI runs, a nonzero fixed seed is expanded by rank as
`base_seed + 1000003 * mpi_rank`. Computational cutoffs are tracked separately
from physical right-censoring so that final evaporation-time files contain only
complete valid unbinding events.

A trajectory can become physically bound only at a scattering. If an
uncaptured trajectory acquires negative energy during scatter-free propagation,
the run classifies that trajectory as a numerical failure instead of allowing
repeated bound Kepler returns to stall its MPI batch.

## Outputs

For non-capture parameter-point runs, the final files are written after MPI
reduction:

- `bincount.txt`: capture-conditioned residence-time and `v^2 dt` radial
  histograms with error estimates. The grid is uniform at `0.001 R_sun`
  through `1.1 R_sun`. From there, shell widths start continuously at
  `0.001 R_sun` and grow geometrically by 2% per shell toward a global
  `10 R_sun` width cap. The default 1-AU domain ends before that cap is reached
  and uses 423 exterior shells. Negative-energy
  exterior Kepler arcs contribute exact shell integrals. A captured orbit whose
  apoapsis exceeds 1 AU is propagated analytically to its outward 1-AU
  crossing, marked `outer_domain_removal`, and contributes capture plus
  residence statistics through that crossing. It does not contribute an
  evaporation-time event. Numerical failures are excluded from residence
  statistics. Accepted numerical RK intervals are conservatively split using
  adaptive Hermite dense output. Uncaptured residence histograms are not
  written because capture membership gates every residence contribution.
- `evaporation_times.txt`: compact complete-event table with
  `rank trajectory_id lifetime_unbinding_sec r_capture_Rsun E_capture_eV
  dE_capture_eV`, followed by the number of negative-energy exterior arcs,
  the first/last/maximum osculating Kepler periods at outward `1.1 R_sun`
  crossings, and the corresponding first/last/maximum analytic exterior
  return times. It is sorted by
  `lifetime_unbinding_sec` with `rank trajectory_id` tie-breakers.
- `residence_jackknife_blocks.tsv`: exactly 64 deterministic blocks assigned by
  `splitmix64(base_seed, rank, trajectory_id) % 64`. Each block contains
  attempted, captured, completed uncaptured escape, accepted residence,
  invalid, and outer-domain-removal counts plus its full radial `dt` and
  `v^2 dt` histograms. The writer refuses to publish the file unless every
  scalar count and every radial bin closes against `bincount.txt`. It is the
  required input for delete-one-block propagation of capture-rate/residence
  covariance into shell and channel flux uncertainties.
- `invalid_trajectories.tsv`: always-on, replayable ledger for trajectories
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

Optional trajectory diagnostics are enabled with `trajectory_summary_enabled = true`.
This adds `run_metadata.json`, `trajectory_summary.tsv`, and
`trajectory_events.tsv` without changing the capture or evaporation state
definitions. Set `trajectory_events_enabled = true` and
`trajectory_trace_rate` in `[0, 1]` to select lifecycle traces by a stable hash;
set `trajectory_trace_seed` to keep that selection identical across runs. The
selection does not consume the physics RNG. Traced trajectories include the
complete pre-initial-condition `std::mt19937` state and shifted initial
condition, plus real-time scatter, state-transition, solar-crossing, escape,
censoring, and numerical-failure events.

The bound-exit period is the point-mass osculating Kepler period inferred from
the negative-energy state at the outward `1.1 R_sun` matching surface. The
exterior elapsed time is the physically used analytic travel time: through
apoapsis to the inbound matching surface for contained arcs, or one-way to the
outward 1-AU removal point for outer-domain arcs. These are kept separate
because the osculating full period includes a point-mass continuation through
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
  provisionally. If its bound exterior arc crosses the 1-AU removal surface,
  the forced publication includes the one-way residence integral through that
  crossing. Status rows remain comments so data readers see only bincount bins.
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
