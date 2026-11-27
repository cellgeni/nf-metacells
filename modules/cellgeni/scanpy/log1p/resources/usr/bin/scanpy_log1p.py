#!/usr/bin/env python3
"""Log-transform a normalised layer, log(1 + x), into layers/log1p.

scanpy log1p on CPU, rapids-singlecell log1p with --gpu. The input layer (--layer-key, default layers/normalized) is
read from --normalized (a step output) or the base object.
"""

import argparse
import sys

import scio

STEP = "log1p"


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--normalized", help="step output with the normalised layer (default: the base object)")
    parser.add_argument("--layer-key", default="layers/normalized", help="element to transform")
    parser.add_argument("--base", type=float, help="logarithm base (default: natural)")
    args = parser.parse_args(argv)
    scio.setup_logging()

    import anndata as ad
    import numpy as np

    reader = scio.Reader(args.adata, [args.normalized])
    adata = ad.AnnData(X=reader.read(args.layer_key).astype(np.float32))
    if args.gpu:
        scio.gpu_init()
        import rapids_singlecell as rsc

        adata.X = scio.to_gpu(adata.X)
        rsc.pp.log1p(adata, base=args.base)
    else:
        import scanpy as sc

        sc.pp.log1p(adata, base=args.base)

    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("layers/log1p", adata.X)
    writer.write(
        params={"layer_key": args.layer_key, "base": args.base, "gpu": args.gpu},
        versions=scio.versions("anndata", "scanpy", "rapids-singlecell-cu12"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
