#!/usr/bin/env python3
"""Normalise counts to the same total per cell into layers/normalized.

scanpy normalize_total on CPU, rapids-singlecell normalize_total with --gpu. Counts are read from --counts-key
(X, layers/<name> or raw/X) of the base object and cast to float32.
"""

import argparse
import sys

import scio

STEP = "normalizetotal"


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--counts-key", default="X", help="element with the counts: X, layers/<name> or raw/X")
    parser.add_argument("--target-sum", type=float, help="total per cell after normalisation (default: median total)")
    parser.add_argument(
        "--exclude-highly-expressed", action="store_true",
        help="leave genes that make up more than --max-fraction of a cell's counts out of its total",
    )
    parser.add_argument("--max-fraction", type=float, default=0.05)
    args = parser.parse_args(argv)
    scio.setup_logging()

    import anndata as ad
    import numpy as np

    reader = scio.Reader(args.adata)
    adata = ad.AnnData(X=reader.read(args.counts_key).astype(np.float32))
    options = {
        "target_sum": args.target_sum,
        "exclude_highly_expressed": args.exclude_highly_expressed,
        "max_fraction": args.max_fraction,
    }
    if args.gpu:
        scio.gpu_init()
        import rapids_singlecell as rsc

        adata.X = scio.to_gpu(adata.X)
        rsc.pp.normalize_total(adata, **options)
    else:
        import scanpy as sc

        sc.pp.normalize_total(adata, **options)

    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("layers/normalized", adata.X)
    writer.write(
        params={"counts_key": args.counts_key, "gpu": args.gpu, **options},
        versions=scio.versions("anndata", "scanpy", "rapids-singlecell-cu12"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
