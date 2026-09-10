# Scripts Overview

This directory contains command-line runners that call high-level workflows from `plot_theory_copy.py`.

## Design goal

- `plot_theory_copy.py` acts as the library / physics backend
- `scripts/run_*.py` act as CLI entry points
- common CLI helpers live in `scripts/_common.py`

---

## Presets

Some runners support:

- `--preset fast`
- `--preset medium`
- `--preset production`

Typical meanings:

- `fast`: quick testing / rough plots
- `medium`: balanced default
- `production`: higher-accuracy runs

You can still override individual parameters manually, e.g.

```bash
--preset fast --n-u 32
