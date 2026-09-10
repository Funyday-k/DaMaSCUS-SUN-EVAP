#!/usr/bin/env python3
import os
import argparse
from _common import (
    BASE_DIR, add_mass_args, add_earth_file_arg, add_output_root_arg,
    add_common_integration_args, parse_masses, load_earth_data,
    call_supported, print_header, print_result_summary
)
import plot_theory_copy as pt


def build_parser():
    parser = argparse.ArgumentParser(description="Run verified-only SD thermal convergence.")
    add_mass_args(parser)
    add_earth_file_arg(parser)
    add_common_integration_args(parser, include_sigmas=("sd",))
    add_output_root_arg(parser, os.path.join(BASE_DIR, "outputs", "sd_thermal_convergence"))
    return parser


def main():
    args = build_parser().parse_args()
    DM_masses = parse_masses(args)
    earth_data = load_earth_data(pt, earth_file=args.earth_file, sd_mode="verified_only")

    print_header(
        "Verified-only SD thermal convergence runner",
        output_root=args.output_root,
        DM_masses=DM_masses,
        sigma_SD_p=args.sigma_sd_p,
    )

    results = call_supported(
        pt.test_verified_only_sd_thermal_convergence,
        earth_data=earth_data,
        DM_masses=DM_masses,
        masses=DM_masses,
        selected_masses=DM_masses,
        sigma_SD_p=args.sigma_sd_p,
        output_root=args.output_root,
        u_max=args.u_max,
        n_u=args.n_u,
        n_t_speed=args.n_t_speed,
        n_t_mu=args.n_t_mu,
        n_scatter_mu=args.n_scatter_mu,
        n_scatter_phi=args.n_scatter_phi,
        max_workers=args.max_workers,
    )
    print_result_summary(results)


if __name__ == "__main__":
    main()
