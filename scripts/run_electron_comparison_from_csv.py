#!/usr/bin/env python3
import os
import argparse
from _common import BASE_DIR, add_output_root_arg, call_supported, print_header, print_result_summary
import plot_theory_copy as pt


def build_parser():
    parser = argparse.ArgumentParser(description="Plot electron comparison from CSV.")
    parser.add_argument("--input-csv", type=str, required=True, help="Input CSV path.")
    add_output_root_arg(parser, os.path.join(BASE_DIR, "outputs", "electron_comparison_from_csv"))
    return parser


def main():
    args = build_parser().parse_args()

    print_header(
        "Electron comparison-from-CSV runner",
        input_csv=args.input_csv,
        output_root=args.output_root,
    )

    results = call_supported(
        pt.plot_electron_comparison_from_csv,
        input_csv=args.input_csv,
        csv_path=args.input_csv,
        source_csv=args.input_csv,
        results_csv=args.input_csv,
        output_root=args.output_root,
    )
    print_result_summary(results)


if __name__ == "__main__":
    main()
