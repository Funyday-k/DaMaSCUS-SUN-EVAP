# DaMaSCUS-SUN-EVAP

Dark Matter Simulation Code for the Sun — Evaporation Branch

## Overview

This project extends [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN) to study dark matter (DM) capture and evaporation in the Sun. It simulates DM particle trajectories inside the Sun using adaptive RK45 integration, tracks gravitational capture via energy monitoring, and computes time-weighted radial distribution histograms (bincount) on-the-fly.

## Key Features

- **Online bincount accumulation**: Replaces trajectory file output with in-memory histogram accumulation (r-histogram and v²-histogram, 2000 bins from 0 to 2R☉)
- **Capture detection**: Monitors particle energy at every integration step; particles with any negative energy are classified as captured
- **Evaporation time**: Records the duration of gravitational binding for captured particles, with truncation flags for censored data
- **MPI parallelization**: Each rank independently targets `ceil(sample_size / N_ranks)` captured particles
- **Safety valve**: `max_trajectories` parameter prevents infinite runtime for low capture-rate scenarios

## Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The executable is placed in `bin/DaMaSCUS-SUN`.

### Dependencies

- C++17 compiler
- MPI (OpenMPI or MPICH)
- CMake ≥ 3.14
- libconfig++

External dependencies (`libphysica`, `obscura`) are fetched automatically by CMake.

## Usage

```bash
mpirun -np 8 bin/DaMaSCUS-SUN bin/config.cfg
```

### Configuration File

Key parameters in the `.cfg` file:

| Parameter | Description |
|-----------|-------------|
| `sample_size` | Target number of captured particles |
| `max_trajectories` | Maximum total trajectories before forced stop (default: sample_size × 1000) |
| `output_dir` | Directory for output files |
| `DM_mass` | Dark matter mass in GeV |
| `DM_cross_section_nucleon` | DM-nucleon cross section in cm² |
| `run_mode` | "Parameter point" or "Parameter scan" |

### Output Files

For each parameter point, the main generated files are:

- `bincount.txt` — Combined captured/not-captured time-weighted radial histogram with error estimates
- `evaporation_summary.txt` — Evaporation time for each captured trajectory, keyed by `rank` and rank-local `trajectory_id`
- `computation_time_summary.txt` — Wall-clock time and RK45 step-count statistics

When `snapshot_enabled = true`, intermediate files are written under `snapshot/`:

- `snapshot_{time}s.txt` — Cumulative snapshot report and bincount histogram. The rank diagnostic table combines the checkpoint source with the current `running`/`done` status.
- `snapshot_{time}s_evaporation.txt` — Rank-aware evaporation list for captured trajectories completed by that snapshot time.

## References

- T. Emken, C. Kouvaris, [DaMaSCUS-SUN](https://github.com/temken/DaMaSCUS-SUN)
- Garani & Palomares-Ruiz (2017), Evaporation of dark matter in the Sun
