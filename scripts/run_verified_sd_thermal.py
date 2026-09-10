#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)

if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

import plot_theory_copy as pt


def parse_masses(args):
    if args.masses is not None:
        masses = [
            float(x.strip())
            for x in args.masses.split(",")
            if x.strip()
        ]
        if len(masses) == 0:
            raise ValueError("--masses was provided but no valid numbers were found.")
        return np.array(masses, dtype=float)

    if args.m_min <= 0.0 or args.m_max <= 0.0:
        raise ValueError("m_min and m_max must be positive.")
    if args.m_max < args.m_min:
        raise ValueError("m_max must be >= m_min.")
    if args.n_masses < 1:
        raise ValueError("n_masses must be >= 1.")

    if args.n_masses == 1:
        return np.array([args.m_min], dtype=float)

    return np.logspace(
        np.log10(args.m_min),
        np.log10(args.m_max),
        args.n_masses
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run verified-only refined SD thermal comparison."
    )

    parser.add_argument(
        "--masses",
        type=str,
        default=None,
        help='Explicit comma-separated DM masses, e.g. "1,10,100".'
    )
    parser.add_argument(
        "--m-min",
        type=float,
        default=0.5,
        help="Minimum DM mass in GeV for logspace grid."
    )
    parser.add_argument(
        "--m-max",
        type=float,
        default=1e3,
        help="Maximum DM mass in GeV for logspace grid."
    )
    parser.add_argument(
        "--n-masses",
        type=int,
        default=40,
        help="Number of DM mass points for logspace grid."
    )

    parser.add_argument(
        "--earth-file",
        type=str,
        default=None,
        help="Optional path to earth_prem.dat. Default: library default."
    )

    parser.add_argument(
        "--sigma-sd-p",
        type=float,
        default=1e-40,
        help="Reference SD proton cross section in cm^2."
    )
    parser.add_argument(
        "--u-max",
        type=float,
        default=800.0,
        help="Maximum halo speed in km/s."
    )
    parser.add_argument(
        "--n-u",
        type=int,
        default=80,
        help="Number of u-grid points."
    )
    parser.add_argument(
        "--n-t-speed",
        type=int,
        default=4,
        help="Thermal speed grid size."
    )
    parser.add_argument(
        "--n-t-mu",
        type=int,
        default=4,
        help="Thermal angular grid size."
    )
    parser.add_argument(
        "--n-scatter-mu",
        type=int,
        default=8,
        help="Scatter mu grid size."
    )
    parser.add_argument(
        "--n-scatter-phi",
        type=int,
        default=12,
        help="Scatter phi grid size."
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Number of worker processes."
    )

    parser.add_argument(
        "--output-root",
        type=str,
        default=os.path.join(BASE_DIR, "outputs", "verified_sd_thermal"),
        help="Root output directory."
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    DM_masses = parse_masses(args)

    print("=" * 90)
    print("Verified-only refined SD thermal runner")
    print("=" * 90)
    print(f"BASE_DIR       = {BASE_DIR}")
    print(f"earth_file     = {args.earth_file}")
    print(f"sd_mode        = verified_only")
    print(f"output_root    = {args.output_root}")
    print(f"DM_masses      = {DM_masses}")
    print(f"sigma_SD_p     = {args.sigma_sd_p:.6e} cm^2")
    print(f"n_u            = {args.n_u}")
    print(f"n_t_speed      = {args.n_t_speed}")
    print(f"n_t_mu         = {args.n_t_mu}")
    print(f"n_scatter_mu   = {args.n_scatter_mu}")
    print(f"n_scatter_phi  = {args.n_scatter_phi}")
    print(f"max_workers    = {args.max_workers}")

    load_kwargs = {"sd_mode": "verified_only"}
    if args.earth_file is not None:
        load_kwargs["filepath"] = args.earth_file

    earth_data = pt.load_earth_composition(**load_kwargs)

    results = pt.plot_verified_only_sd_thermal_constant_refined(
        earth_data=earth_data,
        DM_masses=DM_masses,
        sigma_SD_p=args.sigma_sd_p,
        output_root=args.output_root,
        u_max=args.u_max,
        n_u=args.n_u,
        n_t_speed=args.n_t_speed,
        n_t_mu=args.n_t_mu,
        n_scatter_mu=args.n_scatter_mu,
        n_scatter_phi=args.n_scatter_phi,
        max_workers=args.max_workers
    )

    print("=" * 90)
    print("Run completed")
    print("=" * 90)

    if isinstance(results, dict):
        for key in ("csv_path", "fig_png", "fig_pdf", "log_path"):
            if key in results:
                print(f"{key:10s}: {results[key]}")
        print(f"result keys: {list(results.keys())}")


if __name__ == "__main__":
    main()
