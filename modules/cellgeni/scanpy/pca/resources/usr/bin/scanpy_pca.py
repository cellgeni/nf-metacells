#!/usr/bin/env python3
"""Principal components of a log-normalised layer into obsm/X_pca, varm/PCs and uns/pca.

scanpy pca on CPU, rapids-singlecell pca with --gpu. The layer (--matrix-key, default layers/log1p) comes from
--matrix or the base object and the gene mask (--mask-key, default var/highly_variable) from --mask or the base
object; --no-mask uses every gene. --scale scales the masked genes to unit variance in memory first (as
sc.pp.scale, clipped at --scale-max-value), so no dense layer is written. Loadings of unmasked genes are zero.
"""

import argparse
import sys

import scio

STEP = "pca"


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--matrix", help="step output with the log-normalised layer (default: the base object)")
    parser.add_argument("--matrix-key", default="layers/log1p")
    parser.add_argument("--mask", help="step output with the gene mask (default: the base object)")
    parser.add_argument("--mask-key", default="var/highly_variable")
    parser.add_argument("--no-mask", action="store_true", help="use every gene")
    parser.add_argument("--n-comps", type=int, default=50)
    parser.add_argument("--no-zero-center", action="store_true", help="truncated SVD of the uncentred matrix")
    parser.add_argument("--svd-solver", help="solver name as scanpy/rapids-singlecell know it (default: theirs)")
    parser.add_argument("--scale", action="store_true", help="scale the masked genes to unit variance first")
    parser.add_argument("--scale-max-value", type=float, help="clip scaled values (with --scale)")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    scio.setup_logging()

    import anndata as ad
    import numpy as np

    reader = scio.Reader(args.adata, [args.matrix, args.mask])
    matrix = reader.read(args.matrix_key)
    mask = np.ones(reader.n_vars, dtype=bool) if args.no_mask else reader.read(args.mask_key).to_numpy(dtype=bool)
    n_genes = int(mask.sum())
    n_comps = min(args.n_comps, reader.n_obs - 1, n_genes - 1)
    if n_comps < 1:
        raise SystemExit(f"cannot compute PCA of {reader.n_obs} cells x {n_genes} genes")
    if n_comps < args.n_comps:
        scio.log.warning("%d components requested, %d possible: computing %d", args.n_comps, n_comps, n_comps)

    # Only the masked genes enter the PCA; loadings go back into an n_vars x n_comps matrix
    adata = ad.AnnData(X=matrix[:, np.flatnonzero(mask)].astype(np.float32))
    options = {"n_comps": n_comps, "zero_center": not args.no_zero_center, "svd_solver": args.svd_solver}
    if args.gpu:
        scio.gpu_init()
        import rapids_singlecell as rsc

        adata.X = scio.to_gpu(adata.X)
        if args.scale:
            rsc.pp.scale(adata, zero_center=True, max_value=args.scale_max_value)
        rsc.pp.pca(adata, mask_var=None, rng=args.seed, **options)
    else:
        import scanpy as sc

        if args.scale:
            sc.pp.scale(adata, zero_center=True, max_value=args.scale_max_value)
        sc.pp.pca(adata, mask_var=None, random_state=args.seed, **options)

    loadings = np.zeros((reader.n_vars, n_comps), dtype=np.float32)
    loadings[mask] = scio.to_cpu(adata.varm["PCs"])
    info = adata.uns["pca"]
    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("obsm/X_pca", scio.to_cpu(adata.obsm["X_pca"]).astype(np.float32))
    writer.add("varm/PCs", loadings)
    writer.add("uns/pca", {
        "variance": scio.to_cpu(info["variance"]),
        "variance_ratio": scio.to_cpu(info["variance_ratio"]),
        "params": {"zero_center": options["zero_center"], "mask_var": None if args.no_mask else args.mask_key},
    })
    writer.write(
        params={**options, "matrix_key": args.matrix_key, "mask_key": None if args.no_mask else args.mask_key,
                "n_genes": n_genes, "scale": args.scale, "scale_max_value": args.scale_max_value,
                "seed": args.seed, "gpu": args.gpu},
        versions=scio.versions("anndata", "scanpy", "rapids-singlecell-cu12"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
