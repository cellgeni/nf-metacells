#!/usr/bin/env python3
"""Hierarchical metacells built within groups: `<prefix>.membership.csv` (barcode,metacell).

For each group of --constrain-by, a kNN graph of the embedding (scanpy on CPU, rapids-singlecell with --gpu, which
searches exactly) is cut by a scikit-network hierarchy (Louvain or Paris) into balanced clusters of at most
--n-max cells. A group with fewer than --n-max cells becomes one metacell. The cut is balanced, so clusters below
--n-min are a property of that group's graph rather than something to merge away: they are counted in
`<prefix>.hmetacells.json`, not changed. Cells whose group is missing are left unassigned (an empty metacell).
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd

import scio

METHOD = "hmetacells"


def neighbor_graph(embedding: np.ndarray, n_neighbors: int, seed: int, use_gpu: bool):
    """obsp/connectivities for one group, from the same backends as cellgeni/scanpy/neighbors."""
    import anndata as ad

    adata = ad.AnnData(obs={"_": np.zeros(embedding.shape[0])}, obsm={"X_rep": embedding})
    options = {"n_neighbors": n_neighbors, "use_rep": "X_rep"}
    if use_gpu:
        import rapids_singlecell as rsc

        rsc.pp.neighbors(adata, rng=seed, **options)
    else:
        import scanpy as sc

        sc.pp.neighbors(adata, random_state=seed, **options)
    return adata.obsp["connectivities"]


def split_group(embedding: np.ndarray, n_max: int, method: str, n_neighbors: int, seed: int,
                use_gpu: bool) -> np.ndarray:
    """Cluster indices for one group, balanced to at most n_max cells each."""
    from sknetwork.hierarchy import LouvainHierarchy, Paris, cut_balanced

    # A group can be smaller than the requested neighbourhood; scanpy and rapids-singlecell both refuse that
    k = min(n_neighbors, embedding.shape[0] - 1)
    connectivities = neighbor_graph(embedding, k, seed, use_gpu)
    model = LouvainHierarchy() if method == "louvain" else Paris()
    return cut_balanced(model.fit_transform(connectivities), max_cluster_size=n_max)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--adata", required=True, help="base AnnData object (.h5ad or .zarr)")
    parser.add_argument("--prefix", required=True, help="output name prefix, usually the record id")
    parser.add_argument("--embedding", help="step output with the embedding (default: the base object)")
    parser.add_argument("--embedding-key", default="obsm/X_pca")
    parser.add_argument("--constrain-by", required=True, help="obs column whose groups metacells are built within")
    parser.add_argument("--n-min", type=int, required=True, help="metacell size to report below, not to enforce")
    parser.add_argument("--n-max", type=int, required=True, help="largest metacell the cut may leave")
    parser.add_argument("--method", default="louvain", choices=["louvain", "paris"], help="sknetwork hierarchy")
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--n-dims", type=int, help="use the first n columns of the embedding (default: all)")
    parser.add_argument("--gpu", action="store_true", help="build the kNN graphs on the GPU")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    scio.setup_logging()
    if args.n_min >= args.n_max:
        raise SystemExit(f"--n-min {args.n_min} must be below --n-max {args.n_max}")

    reader = scio.Reader(args.adata, [args.embedding])
    embedding = np.asarray(reader.read(args.embedding_key), dtype=np.float32)
    if args.n_dims:
        embedding = embedding[:, : args.n_dims]
    groups = pd.Series(reader.read(f"obs/{args.constrain_by}")).astype("object")
    if len(groups) != reader.n_obs:
        raise SystemExit(f"obs/{args.constrain_by} has {len(groups)} values for {reader.n_obs} cells")
    if args.gpu:
        scio.gpu_init()

    metacells = np.full(reader.n_obs, "", dtype=object)
    known = groups.notna().to_numpy()
    for group in pd.unique(groups[known]):
        index = np.flatnonzero(known & (groups == group).to_numpy())
        if index.size < args.n_max:
            metacells[index] = f"{group}_0"
            scio.log.info("group %s: %s cells, one metacell", group, index.size)
            continue
        clusters = split_group(embedding[index], args.n_max, args.method, args.n_neighbors, args.seed, args.gpu)
        metacells[index] = [f"{group}_{cluster}" for cluster in clusters]
        scio.log.info("group %s: %s cells into %s metacells", group, index.size, len(set(clusters)))

    membership = pd.DataFrame({"barcode": reader.obs_names, "metacell": metacells})
    membership.to_csv(f"{args.prefix}.membership.csv", index=False)

    sizes = membership.loc[membership["metacell"] != "", "metacell"].value_counts()
    unassigned = int((~known).sum())
    summary = {
        "method": METHOD,
        "n_cells": int(len(membership)),
        "n_metacells": int(sizes.size),
        "n_groups": int(pd.unique(groups[known]).size),
        "unassigned": unassigned,
        "metacell_size": {"min": int(sizes.min()), "median": float(sizes.median()), "max": int(sizes.max())}
        if sizes.size else None,
        # The balanced cut bounds metacells above, not below, so these are reported rather than merged away
        "below_n_min": int((sizes < args.n_min).sum()),
        "constrain_by": args.constrain_by,
        "embedding_key": args.embedding_key,
        "hierarchy": args.method,
        "n_neighbors": args.n_neighbors,
        "n_min": args.n_min,
        "n_max": args.n_max,
        "gpu": bool(args.gpu),
        "seed": args.seed,
        "versions": scio.versions("anndata", "scanpy", "scikit-network", "rapids-singlecell-cu12"),
    }
    with open(f"{args.prefix}.hmetacells.json", "w") as out:
        json.dump(summary, out, indent=2, sort_keys=True)

    if unassigned:
        scio.log.warning("%s cells have no %s and are unassigned", unassigned, args.constrain_by)
    if summary["below_n_min"]:
        scio.log.warning("%s of %s metacells hold fewer than %s cells", summary["below_n_min"], sizes.size,
                         args.n_min)
    scio.log.info("%s metacells for %s cells in %s groups on the %s", sizes.size, len(membership),
                  summary["n_groups"], "GPU" if args.gpu else "CPU")
    return 0


if __name__ == "__main__":
    sys.exit(main())
