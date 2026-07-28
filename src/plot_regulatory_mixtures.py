#!/usr/bin/env python3
"""Plot empirical H3K27ac distributions and guarded Gaussian mixtures."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.regulatory_mixture_plot import build_mixture_distribution_plot


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--mixtures", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    parser.add_argument("--output-bins", type=Path, required=True)
    parser.add_argument("--output-metrics", type=Path, required=True)
    parser.add_argument("--bins", type=int, default=70)
    args = parser.parse_args()
    metrics = build_mixture_distribution_plot(
        catalog_path=args.catalog,
        mixture_path=args.mixtures,
        output_svg=args.output_svg,
        output_bins=args.output_bins,
        output_metrics=args.output_metrics,
        bin_n=args.bins,
    )
    print(
        f"contexts={len(metrics['contexts'])} "
        f"supported={len(metrics['supported_contexts'])} "
        f"unsupported={','.join(metrics['unsupported_contexts']) or 'none'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
