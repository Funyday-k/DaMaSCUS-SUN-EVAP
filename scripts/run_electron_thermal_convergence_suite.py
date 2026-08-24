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
    parser = argparse.ArgumentParser(description="Run electron thermal convergence suite.")
    add_mass_args(parser)
    add_earth_file_arg(parser)
    add_common_integration_args(parser, include_sigmas=("electron",))
    add_output_root_arg(parser, os.path.join(BASE_DIR, "outputs", "electron_thermal_convergence_suite"))
    return parser


def main():
    args = build_parser().parse_args()
    DM_masses = parse_masses(args)
    earth_data = load_earth_data(pt, earth_file=args.earth_file)

    print_header(
        "Electron thermal convergence suite runner",
        output_root=args.output_root,
        DM_masses=DM_masses,
        sigma_electron=args.sigma_electron,
    )

    results = call_supported(
        pt.run_and_plot_electron_thermal_convergence_suite,
        earth_data=earth_data,
        DM_masses=DM_masses,
        sigma_electron=args.sigma_electron,
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
