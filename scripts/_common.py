#!/usr/bin/env python3
import os
import sys
import inspect
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)

if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)


def add_mass_args(parser, m_min=0.5, m_max=1e3, n_masses=40):
    parser.add_argument(
        "--masses",
        type=str,
        default=None,
        help='Explicit comma-separated DM masses, e.g. "1,10,100".',
    )
    parser.add_argument(
        "--m-min",
        type=float,
        default=m_min,
        help="Minimum DM mass in GeV for logspace grid.",
    )
    parser.add_argument(
        "--m-max",
        type=float,
        default=m_max,
        help="Maximum DM mass in GeV for logspace grid.",
    )
    parser.add_argument(
        "--n-masses",
        type=int,
        default=n_masses,
        help="Number of DM mass points for logspace grid.",
    )


def add_earth_file_arg(parser):
    parser.add_argument(
        "--earth-file",
        type=str,
        default=None,
        help="Optional path to earth_prem.dat. Default: library default.",
    )


def add_output_root_arg(parser, default_root):
    parser.add_argument(
        "--output-root",
        type=str,
        default=default_root,
        help="Root output directory.",
    )


def add_common_integration_args(
    parser,
    include_sigmas=("sd",),
    include_shell_step=False,
    include_u_grid_mode=False,
    defaults=None,
):
    if defaults is None:
        defaults = {}

    if "si" in include_sigmas:
        parser.add_argument(
            "--sigma-si-p",
            type=float,
            default=defaults.get("sigma_si_p", 1e-40),
            help="Reference SI proton cross section in cm^2.",
        )

    if "sd" in include_sigmas:
        parser.add_argument(
            "--sigma-sd-p",
            type=float,
            default=defaults.get("sigma_sd_p", 1e-40),
            help="Reference SD proton cross section in cm^2.",
        )

    if "electron" in include_sigmas:
        parser.add_argument(
            "--sigma-electron",
            type=float,
            default=defaults.get("sigma_electron", 1e-40),
            help="Reference electron cross section in cm^2.",
        )

    parser.add_argument(
        "--u-max",
        type=float,
        default=defaults.get("u_max", 800.0),
        help="Maximum halo speed in km/s.",
    )
    parser.add_argument(
        "--n-u",
        type=int,
        default=defaults.get("n_u", 80),
        help="Number of u-grid points.",
    )
    parser.add_argument(
        "--n-t-speed",
        type=int,
        default=defaults.get("n_t_speed", 4),
        help="Thermal speed grid size.",
    )
    parser.add_argument(
        "--n-t-mu",
        type=int,
        default=defaults.get("n_t_mu", 4),
        help="Thermal angular grid size.",
    )
    parser.add_argument(
        "--n-scatter-mu",
        type=int,
        default=defaults.get("n_scatter_mu", 8),
        help="Scatter mu grid size.",
    )
    parser.add_argument(
        "--n-scatter-phi",
        type=int,
        default=defaults.get("n_scatter_phi", 12),
        help="Scatter phi grid size.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=defaults.get("max_workers", None),
        help="Number of worker processes.",
    )

    if include_shell_step:
        parser.add_argument(
            "--shell-step",
            type=int,
            default=defaults.get("shell_step", 1),
            help="Radial shell stride.",
        )

    if include_u_grid_mode:
        parser.add_argument(
            "--u-grid-mode",
            type=str,
            default=defaults.get("u_grid_mode", "log"),
            choices=["linear", "log", "hybrid"],
            help="Halo speed grid mode.",
        )


def parse_masses(args):
    masses_str = getattr(args, "masses", None)
    if masses_str is not None:
        masses = [float(x.strip()) for x in masses_str.split(",") if x.strip()]
        if not masses:
            raise ValueError("--masses was provided but no valid numbers were found.")
        return np.array(masses, dtype=float)

    m_min = getattr(args, "m_min")
    m_max = getattr(args, "m_max")
    n_masses = getattr(args, "n_masses")

    if m_min <= 0.0 or m_max <= 0.0:
        raise ValueError("m_min and m_max must be positive.")
    if m_max < m_min:
        raise ValueError("m_max must be >= m_min.")
    if n_masses < 1:
        raise ValueError("n_masses must be >= 1.")

    if n_masses == 1:
        return np.array([m_min], dtype=float)

    return np.logspace(np.log10(m_min), np.log10(m_max), n_masses)


def load_earth_data(pt, earth_file=None, sd_mode=None, min_mass_fraction=None):
    kwargs = {}
    if earth_file is not None:
        kwargs["filepath"] = earth_file
    if sd_mode is not None:
        kwargs["sd_mode"] = sd_mode
    if min_mass_fraction is not None:
        kwargs["min_mass_fraction"] = min_mass_fraction
    return pt.load_earth_composition(**kwargs)


def call_supported(func, **kwargs):
    sig = inspect.signature(func)
    supported = {k: v for k, v in kwargs.items() if k in sig.parameters}
    return func(**supported)


def print_header(title, **items):
    line = "=" * 90
    print(line)
    print(title)
    print(line)
    for k, v in items.items():
        print(f"{k:16s}= {v}")


def print_result_summary(results):
    line = "=" * 90
    print(line)
    print("Run completed")
    print(line)

    if isinstance(results, dict):
        for key in (
            "csv_path",
            "fig_png",
            "fig_pdf",
            "log_path",
            "summary_csv",
            "output_csv",
            "results_csv",
        ):
            if key in results:
                print(f"{key:16s}= {results[key]}")
        print(f"result keys      = {list(results.keys())}")
