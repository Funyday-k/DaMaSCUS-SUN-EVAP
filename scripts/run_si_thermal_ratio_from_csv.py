#!/usr/bin/env python3
import os
import argparse
from _common import BASE_DIR, add_output_root_arg, call_supported, print_header, print_result_summary
import plot_theory_copy as pt


def build_parser():
    parser = argparse.ArgumentParser(description="Plot SI thermal-ratio comparison from CSV.")
    parser.add_argument("--thermal-csv", type=str, required=True, help="Thermal CSV path.")
    parser.add_argument("--constant-csv", type=str, required=True, help="Constant/T0 CSV path.")
    add_output_root_arg(parser, os.path.join(BASE_DIR, "outputs", "si_thermal_ratio_from_csv"))
    return parser


def main():
    args = build_parser().parse_args()

    print_header(
        "SI thermal-ratio-from-CSV runner",
        thermal_csv=args.thermal_csv,
        constant_csv=args.constant_csv,
        output_root=args.output_root,
    )

    results = call_supported(
        pt.plot_si_thermal_ratio_comparison_from_csv,
        thermal_csv=args.thermal_csv,
        constant_csv=args.constant_csv,
        t0_csv=args.constant_csv,
        baseline_csv=args.constant_csv,
        output_root=args.output_root,
    )
    print_result_summary(results)


if __name__ == "__main__":
    main()
