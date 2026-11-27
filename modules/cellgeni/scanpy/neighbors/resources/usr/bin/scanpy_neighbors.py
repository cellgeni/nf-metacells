#!/usr/bin/env python3
"""k-nearest-neighbour graph of an embedding into obsp/connectivities, obsp/distances and uns/neighbors.

scanpy neighbors on CPU, rapids-singlecell neighbors with --gpu (exact brute-force search by default; scanpy uses
an approximate search on large data). The embedding (--embedding-key, default obsm/X_pca) comes from --embedding or
the base object.
"""

import argparse
import sys

import scio

STEP = "neighbors"


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--embedding", help="step output with the embedding (default: the base object)")
    parser.add_argument("--embedding-key", default="obsm/X_pca")
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--n-dims", type=int, help="use the first n columns of the embedding (default: all)")
    parser.add_argument("--metric", default="euclidean")
    parser.add_argument("--method", default="umap", choices=["umap", "gauss"], help="connectivity kernel")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    scio.setup_logging()

    import anndata as ad
    import numpy as np

    reader = scio.Reader(args.adata, [args.embedding])
    embedding = np.asarray(reader.read(args.embedding_key), dtype=np.float32)
    if args.n_dims:
        embedding = embedding[:, : args.n_dims]
    if args.n_neighbors >= reader.n_obs:
        raise SystemExit(f"--n-neighbors {args.n_neighbors} needs more than {reader.n_obs} cells")
    adata = ad.AnnData(obs={"_": np.zeros(reader.n_obs)}, obsm={"X_rep": embedding})
    options = {"n_neighbors": args.n_neighbors, "use_rep": "X_rep", "metric": args.metric, "method": args.method}
    if args.gpu:
        scio.gpu_init()
        import rapids_singlecell as rsc

        rsc.pp.neighbors(adata, rng=args.seed, **options)
    else:
        import scanpy as sc

        sc.pp.neighbors(adata, random_state=args.seed, **options)

    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("obsp/connectivities", adata.obsp["connectivities"])
    writer.add("obsp/distances", adata.obsp["distances"])
    writer.add("uns/neighbors", {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": args.n_neighbors, "method": args.method, "metric": args.metric,
                   "use_rep": args.embedding_key.split("/", 1)[1], "random_state": args.seed},
    })
    writer.write(
        params={**options, "embedding_key": args.embedding_key, "n_dims": args.n_dims, "seed": args.seed,
                "gpu": args.gpu},
        versions=scio.versions("anndata", "scanpy", "rapids-singlecell-cu12"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
