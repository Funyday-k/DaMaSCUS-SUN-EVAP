#!/usr/bin/env python3
import os
import argparse
from _common import BASE_DIR, add_output_root_arg, call_supported, print_header, print_result_summary
import plot_theory_copy as pt


def build_parser():
    parser = argparse.ArgumentParser(description="Redraw SI q^2 result from CSV.")
    parser.add_argument("--input-csv", type=str, required=True, help="Input CSV path.")
    add_output_root_arg(parser, os.path.join(BASE_DIR, "outputs", "si_q2_redraw"))
    return parser


def main():
    args = build_parser().parse_args()

    print_header(
        "SI q^2 redraw-from-CSV runner",
        input_csv=args.input_csv,
        output_root=args.output_root,
    )

    results = call_supported(
        pt.redraw_q2_si_from_csv,
        input_csv=args.input_csv,
        csv_path=args.input_csv,
        source_csv=args.input_csv,
        results_csv=args.input_csv,
        output_root=args.output_root,
    )
    print_result_summary(results)


if __name__ == "__main__":
    main()
