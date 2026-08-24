#!/usr/bin/env python3
import os
import argparse
from _common import BASE_DIR, add_output_root_arg, call_supported, print_header, print_result_summary
import plot_theory_copy as pt


def build_parser():
    parser = argparse.ArgumentParser(description="Plot electron thermal convergence summary.")
    parser.add_argument("--summary-csv", type=str, required=True, help="Summary CSV path.")
    add_output_root_arg(parser, os.path.join(BASE_DIR, "outputs", "electron_thermal_convergence_summary"))
    return parser


def main():
    args = build_parser().parse_args()

    print_header(
        "Electron thermal convergence summary runner",
        summary_csv=args.summary_csv,
        output_root=args.output_root,
    )

    results = call_supported(
        pt.plot_electron_thermal_convergence_summary,
        summary_csv=args.summary_csv,
        input_csv=args.summary_csv,
        csv_path=args.summary_csv,
        results_csv=args.summary_csv,
        output_root=args.output_root,
    )
    print_result_summary(results)


if __name__ == "__main__":
    main()
