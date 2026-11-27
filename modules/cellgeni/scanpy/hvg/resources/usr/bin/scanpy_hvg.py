#!/usr/bin/env python3
"""Flag highly variable genes in var/highly_variable, with the flavour's statistics in other var columns.

scanpy highly_variable_genes on CPU, rapids-singlecell with --gpu. Never subsets the genes. Flavours seurat and
cell_ranger read the log-normalised layer (--log1p, --layer-key); seurat_v3 and seurat_v3_paper read the counts of
the base object (--counts-key). --batch-key names an obs column of the base object.
"""

import argparse
import sys

import scio

STEP = "hvg"
COUNT_FLAVORS = ("seurat_v3", "seurat_v3_paper")


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--log1p", help="step output with the log-normalised layer (default: the base object)")
    parser.add_argument("--layer-key", default="layers/log1p", help="log-normalised element (seurat, cell_ranger)")
    parser.add_argument("--counts-key", default="X", help="counts element (seurat_v3 flavours)")
    parser.add_argument("--flavor", default="seurat", choices=["seurat", "cell_ranger", *COUNT_FLAVORS])
    parser.add_argument("--n-top-genes", type=int, help="number of genes to flag (default: cut-offs instead)")
    parser.add_argument("--batch-key", help="obs column: select genes within each batch, then merge")
    parser.add_argument("--min-mean", type=float, default=0.0125)
    parser.add_argument("--max-mean", type=float, default=3)
    parser.add_argument("--min-disp", type=float, default=0.5)
    parser.add_argument("--span", type=float, default=0.3, help="loess span (seurat_v3 flavours)")
    parser.add_argument("--n-bins", type=int, default=20)
    args = parser.parse_args(argv)
    scio.setup_logging()

    import anndata as ad
    import numpy as np

    reader = scio.Reader(args.adata, [args.log1p])
    key = args.counts_key if args.flavor in COUNT_FLAVORS else args.layer_key
    adata = ad.AnnData(X=reader.read(key).astype(np.float32))
    adata.obs_names, adata.var_names = reader.obs_names, reader.var_names
    if args.batch_key:
        adata.obs[args.batch_key] = reader.read(f"obs/{args.batch_key}").to_numpy()
    options = {
        "flavor": args.flavor,
        "n_top_genes": args.n_top_genes,
        "batch_key": args.batch_key,
        "min_mean": args.min_mean,
        "max_mean": args.max_mean,
        "min_disp": args.min_disp,
        "span": args.span,
        "n_bins": args.n_bins,
    }
    if args.gpu:
        scio.gpu_init()
        import rapids_singlecell as rsc

        adata.X = scio.to_gpu(adata.X)
        rsc.pp.highly_variable_genes(adata, **options)
    else:
        import scanpy as sc

        sc.pp.highly_variable_genes(adata, **options)

    selected = int(adata.var["highly_variable"].sum())
    if selected < 2:
        raise SystemExit(f"only {selected} highly variable genes: relax the cut-offs or set --n-top-genes")
    writer = scio.writer_from_args(reader, STEP, args)
    for column in adata.var.columns:
        writer.add(f"var/{column}", adata.var[column])
    writer.add("uns/hvg", {"flavor": args.flavor, "n_selected": selected})
    writer.write(
        params={**options, "input_key": key, "gpu": args.gpu},
        versions=scio.versions("anndata", "scanpy", "rapids-singlecell-cu12"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
