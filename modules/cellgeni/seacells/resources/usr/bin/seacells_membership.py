#!/usr/bin/env python3
"""SEACells metacells from a precomputed embedding: `<prefix>.membership.csv` (barcode,metacell).

Reads obsm/<--embedding-key> through scio (a raw step output, a slim object, or the base object itself), fits the
SEACells archetypal model and writes the hard assignments. The unified model (SEACells a1173c8) runs on CPU and,
with --gpu, on the device; --use-sparse asks for the legacy sparse CPU backend instead, which SEACells refuses to
combine with the unified one. Non-convergence is flagged in `<prefix>.seacells.json`, never raised: SEACells sets
the labels before it complains, and with a seed a rerun reproduces the same partition, so a retry cannot help. A
cuSPARSE failure on the GPU falls back to the CPU once. Randomness comes from the global `np.random`, so --seed
fixes a run on one backend; CPU and GPU give different but equally good partitions at scale, which is why parity
is judged on RSS, purity and metacell sizes rather than ARI.
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd

import scio

METHOD = "seacells"
# cuSPARSE SpGEMM has failed once in the kernel step on a rerun of the same data (CUSPARSE_STATUS_INTERNAL_ERROR)
CUSPARSE_MARKERS = ("CUSPARSE", "cusparse", "cuSPARSE")


def positive_neighbor_counts(embedding: np.ndarray, n_neighbors: int) -> np.ndarray:
    """Non-zero distance counts in scanpy's neighbour graph; scipy drops zero-distance edges."""
    import scanpy as sc

    temp = sc.AnnData(embedding)
    sc.pp.neighbors(temp, n_pcs=0, n_neighbors=n_neighbors)
    return np.diff(temp.obsp["distances"].tocsr().indptr)


def jitter_exact_ties(embedding: np.ndarray, n_neighbors: int = 15, palantir_knn: int = 30,
                      seed: int = 0) -> tuple[np.ndarray, float]:
    """Break exact ties that leave SEACells and Palantir with empty kNN rows (kept from the v1 module).

    Identical coordinates give exactly zero neighbour distances, which sparse matrices do not store: Palantir then
    indexes past the end of a row and SEACells builds a zero kernel bandwidth and collapses into one metacell. The
    tied coordinates carry no information those graph builders can use, so jitter them deterministically, and only
    when a preflight graph shows too few positive distances. Returns the embedding and the jitter scale (0 = none).
    """
    if embedding.ndim != 2:
        raise ValueError(f"the embedding must be 2D, got shape {embedding.shape}")
    if not np.all(np.isfinite(embedding)):
        raise ValueError(f"the embedding has {int((~np.isfinite(embedding)).sum())} non-finite values")
    n_obs = embedding.shape[0]
    if n_obs <= 2:
        return embedding, 0.0
    knn = min(palantir_knn, n_obs - 1)
    if knn < 2:
        return embedding, 0.0
    # Palantir takes adaptive_k = floor(knn / 3); SEACells takes k // 2. Satisfy both.
    required = max(1, min(knn, palantir_knn // 3), min(knn, n_neighbors // 2))
    if not np.any(positive_neighbor_counts(embedding, knn) < required):
        return embedding, 0.0

    deviations = np.std(embedding, axis=0)
    deviations = deviations[np.isfinite(deviations) & (deviations > 0)]
    if deviations.size:
        scale = float(np.median(deviations))
    else:
        finite = embedding[np.isfinite(embedding)]
        scale = float(np.max(np.abs(finite))) if finite.size and np.max(np.abs(finite)) > 0 else 1.0
    jitter = max(scale * 1e-6, float(np.finfo(np.float64).eps))
    jittered = embedding.astype(np.float64, copy=True)
    jittered += np.random.RandomState(seed).normal(0.0, jitter, size=jittered.shape)
    scio.log.warning("jittered %s tied cells by %.3g to keep the kNN graph non-empty", n_obs, jitter)
    return jittered, jitter


def metacell_count(n_obs: int, n_metacells: int | None, gamma: int | None) -> int:
    """--n-metacells wins; otherwise one metacell per --gamma cells, at least one and at most one per cell."""
    if n_metacells is None:
        if gamma is None or gamma <= 0:
            raise SystemExit("give --n-metacells or a positive --gamma")
        n_metacells = max(1, round(n_obs / gamma))
    if not 1 <= n_metacells <= n_obs:
        raise SystemExit(f"--n-metacells {n_metacells} is outside 1..{n_obs} (the number of cells)")
    return int(n_metacells)


def build_model(adata, key: str, n_metacells: int, args: argparse.Namespace, use_gpu: bool):
    """A fitted SEACells model and whether it converged; the labels are set either way."""
    import SEACells

    # Palantir needs 2 <= n_waypoint_eigs <= n_metacells + 1: one eigenvector gives no diffusion components and
    # too many make int(num_waypoints / n_components) zero.
    eigs = min(max(2, min(args.n_waypoint_eigs, n_metacells + 1)), n_metacells + 1)
    if eigs != args.n_waypoint_eigs:
        scio.log.warning("n_waypoint_eigs %s -> %s for %s metacells", args.n_waypoint_eigs, eigs, n_metacells)
    model = SEACells.core.SEACells(
        adata,
        build_kernel_on=key,
        n_SEACells=n_metacells,
        use_gpu=use_gpu,
        use_unified=not args.use_sparse,
        use_sparse=args.use_sparse,
        n_waypoint_eigs=eigs,
        n_neighbors=args.n_neighbors,
        convergence_epsilon=args.convergence_epsilon,
        verbose=False,
    )
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    # Both backends *raise* RuntimeWarning, and only after setting the labels, so catching it as an exception
    # leaves a usable partition. Do not turn warnings into errors here: SEACells and numpy also warn with
    # RuntimeWarning (zero kernel bandwidth, divide by zero), which would be misreported as non-convergence.
    try:
        model.fit(min_iter=args.min_iter, max_iter=args.max_iter)
    except RuntimeWarning as stopped:
        scio.log.warning("SEACells did not converge in %s iterations: %s", args.max_iter, stopped)
        return model, False
    return model, True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--adata", required=True, help="base AnnData object (.h5ad or .zarr)")
    parser.add_argument("--prefix", required=True, help="output name prefix, usually the record id")
    parser.add_argument("--embedding", help="step output holding the embedding (default: the base object)")
    parser.add_argument("--embedding-key", default="X_pca", help="obsm key of the embedding")
    parser.add_argument("--gpu", action="store_true", help="fit on the GPU (unified model)")
    size = parser.add_mutually_exclusive_group(required=True)
    size.add_argument("--n-metacells", type=int, help="number of metacells")
    size.add_argument("--gamma", type=int, help="cells per metacell")
    parser.add_argument("--n-waypoint-eigs", type=int, default=10)
    parser.add_argument("--n-neighbors", type=int, default=15, help="kernel graph neighbours")
    parser.add_argument("--convergence-epsilon", type=float, default=1e-5)
    parser.add_argument("--min-iter", type=int, default=10)
    parser.add_argument("--max-iter", type=int, default=50)
    parser.add_argument("--use-sparse", action="store_true", help="legacy sparse CPU backend (not with --gpu)")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    scio.setup_logging()
    if args.gpu and args.use_sparse:
        raise SystemExit("--use-sparse is a legacy CPU backend: drop it or drop --gpu")

    import anndata as ad

    reader = scio.Reader(args.adata, [args.embedding])
    key = args.embedding_key
    embedding = np.asarray(reader.read(f"obsm/{key}"))
    np.random.seed(args.seed)
    embedding, jitter = jitter_exact_ties(embedding, args.n_neighbors, seed=args.seed)
    n_metacells = metacell_count(len(reader.obs_names), args.n_metacells, args.gamma)

    adata = ad.AnnData(obs=pd.DataFrame(index=pd.Index(reader.obs_names, dtype=str)))
    adata.obsm[key] = embedding
    if args.gpu:
        scio.gpu_init()
    used_gpu = args.gpu
    try:
        model, converged = build_model(adata, key, n_metacells, args, used_gpu)
    except Exception as failed:  # noqa: BLE001 - the message is the only handle cuSPARSE gives us
        if not (args.gpu and any(marker in str(failed) for marker in CUSPARSE_MARKERS)):
            raise
        scio.log.warning("cuSPARSE failed (%s); falling back to the CPU", failed)
        used_gpu = False
        np.random.seed(args.seed)
        model, converged = build_model(adata, key, n_metacells, args, use_gpu=False)

    labels = model.get_hard_assignments()
    membership = pd.DataFrame({"barcode": reader.obs_names, "metacell": labels["SEACell"].to_numpy()})
    membership.to_csv(f"{args.prefix}.membership.csv", index=False)

    sizes = membership["metacell"].value_counts()
    summary = {
        "method": METHOD,
        "converged": converged,
        "n_cells": int(len(membership)),
        "n_metacells": int(sizes.size),
        "requested_metacells": n_metacells,
        "metacell_size": {"min": int(sizes.min()), "median": float(sizes.median()), "max": int(sizes.max())},
        "rss": float(model.compute_RSS()) if hasattr(model, "compute_RSS") else None,
        "embedding_key": key,
        "backend": "sparse" if args.use_sparse else "unified",
        # The device that produced this partition, not the one asked for: a cuSPARSE fallback lands on the CPU,
        # and the two give different partitions, so a parity check must not read a CPU result as a GPU one
        "gpu": bool(used_gpu),
        "gpu_fallback": bool(args.gpu and not used_gpu),
        "jitter": jitter,
        "seed": args.seed,
        "max_iter": args.max_iter,
        "versions": scio.versions("SEACells", "anndata", "scanpy"),
    }
    with open(f"{args.prefix}.seacells.json", "w") as out:
        json.dump(summary, out, indent=2, sort_keys=True)

    scio.log.info("%s metacells for %s cells on the %s (converged=%s)",
                  sizes.size, len(membership), "GPU" if used_gpu else "CPU", converged)
    return 0


if __name__ == "__main__":
    sys.exit(main())
