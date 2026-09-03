#!/usr/bin/env python3
import argparse
from _common import call_supported, print_header
import plot_theory_copy as pt


def build_parser():
    parser = argparse.ArgumentParser(description="Print electron thermal recommended production block.")
    parser.add_argument("--summary-csv", type=str, default=None, help="Optional summary CSV path.")
    return parser


def main():
    args = build_parser().parse_args()

    print_header(
        "Electron thermal recommendation runner",
        summary_csv=args.summary_csv,
    )

    results = call_supported(
        pt.print_electron_thermal_recommended_production_block,
        summary_csv=args.summary_csv,
        input_csv=args.summary_csv,
        csv_path=args.summary_csv,
        results_csv=args.summary_csv,
    )

    if results is not None:
        print(results)


if __name__ == "__main__":
    main()
