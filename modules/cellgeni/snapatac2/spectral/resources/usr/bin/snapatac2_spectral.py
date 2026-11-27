#!/usr/bin/env python3
"""Spectral embedding of peak counts into obsm/X_spectral and uns/spectral_eigenvalue (Laplacian eigenmaps).

CPU: snapATAC2 spectral (matrix-free cosine similarity of IDF-weighted counts; Nystrom with --sample-size; jaccard
with --distance-metric jaccard, which needs --sample-size on large data). GPU (--gpu): a CuPy port of snapATAC2's
matrix-free cosine algorithm (IDF log(n / (1 + df)), row L2 normalisation, degree normalisation, eigsh on
X~ X~^T - D^-1), cosine without sampling only. Features: --features-key mask (default var/selected, from --selected
or the base object) or --all-features. With weighted_by_sd (default) only positive eigenvalues are kept and the
vectors are scaled by their square roots, so the embedding may have fewer than --n-comps columns.
"""

import argparse
import sys

import scio

STEP = "spectral"


ROW_SCALE = r"""
extern "C" __global__ void row_scale(const int n_rows, const int* indptr, float* data, const float* scales) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row < n_rows) {
        const float s = scales[row];
        for (int j = indptr[row]; j < indptr[row + 1]; ++j) data[j] *= s;
    }
}
"""


def spectral_gpu(counts, n_comps, weighted_by_sd=True, random_state=0, eps=1e-8):
    """CuPy port of snapATAC2's cosine spectral embedding (from cellgeni/nf-cluster snapatac2/cupyspectral)."""
    import cupy as cp
    from cupyx.scipy.sparse.linalg import LinearOperator, eigsh

    kernel = cp.RawKernel(ROW_SCALE, "row_scale")

    def scale_rows(m, scales):
        kernel(((m.shape[0] + 255) // 256,), (256,),
               (cp.int32(m.shape[0]), m.indptr.astype(cp.int32), m.data, cp.asarray(scales, dtype=cp.float32)))

    x = scio.to_gpu(counts, dtype=cp.float32).tocsr()
    x.sum_duplicates()
    x.eliminate_zeros()
    n_cells, n_features = x.shape
    # IDF as snapATAC2: log(n / (1 + number of cells with the feature))
    df = cp.bincount(x.indices, minlength=n_features).astype(cp.float32)
    x.data *= cp.log(cp.float32(n_cells) / (1.0 + df))[x.indices]
    row_norm = cp.sqrt(cp.maximum(cp.asarray(x.multiply(x).sum(axis=1)).ravel(), 1e-12))
    scale_rows(x, 1.0 / row_norm)
    # W = Xn Xn^T - I, degree d = W 1; D^-1/2 W D^-1/2 = X~ X~^T - D^-1 with X~ = D^-1/2 Xn
    degree = cp.maximum(cp.asarray(x @ (x.T @ cp.ones(n_cells, dtype=cp.float32))).ravel() - 1.0, eps)
    scale_rows(x, 1.0 / cp.sqrt(degree))
    inv_degree = (1.0 / degree).astype(cp.float32)
    xt = x.T.tocsr()

    def matvec(v):
        v = cp.asarray(v, dtype=cp.float32).ravel()
        return x @ (xt @ v) - inv_degree * v

    operator = LinearOperator((n_cells, n_cells), matvec=matvec, dtype=cp.float32)
    v0 = cp.random.RandomState(random_state).standard_normal(n_cells).astype(cp.float32)
    evals, evecs = eigsh(operator, k=n_comps, which="LM", v0=v0 / (cp.linalg.norm(v0) + eps))
    order = cp.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    if weighted_by_sd:
        keep = evals > 0
        evals = evals[keep]
        evecs = evecs[:, keep] * cp.sqrt(evals)
    return scio.to_cpu(evals), scio.to_cpu(evecs)


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--counts-key", default="X", help="element with the counts: X, layers/<name> or raw/X")
    parser.add_argument("--selected", help="step output with the feature mask (default: the base object)")
    parser.add_argument("--features-key", default="var/selected")
    parser.add_argument("--all-features", action="store_true", help="use every feature")
    parser.add_argument("--n-comps", type=int, default=30)
    parser.add_argument("--distance-metric", default="cosine", choices=["cosine", "jaccard"])
    parser.add_argument("--no-weighted-by-sd", action="store_true", help="keep every eigenvector, unscaled")
    parser.add_argument("--sample-size", type=float, help="Nystrom landmark cells: a count, or a fraction <= 1")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    scio.setup_logging()
    if args.gpu and (args.distance_metric != "cosine" or args.sample_size):
        raise SystemExit("the GPU path computes the full cosine embedding only: drop --sample-size and jaccard")

    import numpy as np

    reader = scio.Reader(args.adata, [args.selected])
    counts = reader.read(args.counts_key)
    mask = np.ones(reader.n_vars, dtype=bool) if args.all_features else reader.read(args.features_key).to_numpy(bool)
    # As snapATAC2: at most n_vars - 1 and n_obs - 1 components
    n_comps = min(args.n_comps, reader.n_vars - 1, reader.n_obs - 1)
    weighted = not args.no_weighted_by_sd
    if args.gpu:
        scio.gpu_init()
        evals, evecs = spectral_gpu(counts[:, np.flatnonzero(mask)], n_comps, weighted, args.seed)
    else:
        import anndata as ad
        import scipy.sparse as sp
        import snapatac2 as snap

        sample_size = None if args.sample_size is None else (
            args.sample_size if args.sample_size <= 1 else int(args.sample_size))
        # snapATAC2's Rust backend takes CSR only; real ATAC files are often stored CSC. Convert before the cast,
        # so the float64 copy is made once and in the layout snapATAC2 wants.
        if sp.issparse(counts) and counts.format != "csr":
            counts = counts.tocsr()
        adata = ad.AnnData(X=counts.astype(np.float64))
        evals, evecs = snap.tl.spectral(
            adata, n_comps=n_comps, features=mask, random_state=args.seed, sample_size=sample_size,
            distance_metric=args.distance_metric, weighted_by_sd=weighted, inplace=False, num_threads=args.threads,
        )
        # The Nystrom path orthogonalises in complex arithmetic; the imaginary parts are zero
        evals, evecs = np.real(evals), np.real(evecs)
    if evecs.shape[1] < n_comps:
        scio.log.warning("%d of %d components kept (non-positive eigenvalues dropped)", evecs.shape[1], n_comps)

    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("obsm/X_spectral", np.asarray(evecs, dtype=np.float32))
    writer.add("uns/spectral_eigenvalue", np.asarray(evals, dtype=np.float64))
    writer.write(
        params={"counts_key": args.counts_key, "features_key": None if args.all_features else args.features_key,
                "n_features": int(mask.sum()), "n_comps": n_comps, "distance_metric": args.distance_metric,
                "weighted_by_sd": weighted, "sample_size": args.sample_size, "seed": args.seed, "gpu": args.gpu},
        versions=scio.versions("anndata", "snapatac2", "cupy-cuda12x"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
