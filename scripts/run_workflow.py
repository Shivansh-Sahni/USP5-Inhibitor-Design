from __future__ import annotations

import argparse
from pathlib import Path

from usp5_workflow.pipeline import run_workflow


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the USP5 inhibitor modeling baseline workflow.")
    parser.add_argument("--input", type=Path, default=Path("data/raw/First.csv"), help="Input CSV file.")
    parser.add_argument("--output", type=Path, default=Path("outputs"), help="Output directory.")
    parser.add_argument(
        "--fingerprint-sizes",
        nargs="+",
        type=int,
        default=[512, 1024, 2048],
        help="Morgan fingerprint bit sizes to evaluate.",
    )
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed for reproducibility.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_workflow(
        input_path=args.input,
        output_dir=args.output,
        fingerprint_sizes=args.fingerprint_sizes,
        random_seed=args.random_seed,
    )


if __name__ == "__main__":
    main()
