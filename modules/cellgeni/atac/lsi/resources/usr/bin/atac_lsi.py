#!/usr/bin/env python3
"""Latent semantic indexing of a TF-IDF layer into obsm/X_lsi, varm/LSI and uns/lsi, as muon's lsi does.

Truncated SVD of --matrix-key (default layers/tfidf, from --tfidf or the base object): scipy svds in float64 on CPU,
cupyx svds with --gpu. Components are ordered by decreasing singular value and signed so that each cell vector's
largest entry is positive (sklearn's svd_flip; muon leaves the sign to the solver). Cell embeddings are z-scored per
component unless --no-scale-embeddings; --drop-first drops the first component, which mostly tracks sequencing
depth (muon keeps it). uns/lsi holds stdev = singular value / sqrt(n_cells - 1).
"""

import argparse
import sys

import scio

STEP = "lsi"


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--tfidf", help="step output with the TF-IDF layer (default: the base object)")
    parser.add_argument("--matrix-key", default="layers/tfidf")
    parser.add_argument("--n-comps", type=int, default=50)
    parser.add_argument("--no-scale-embeddings", action="store_true", help="keep the raw left singular vectors")
    parser.add_argument("--drop-first", action="store_true", help="drop the first component")
    args = parser.parse_args(argv)
    scio.setup_logging()

    import numpy as np

    reader = scio.Reader(args.adata, [args.tfidf])
    matrix = reader.read(args.matrix_key)
    n_comps = min(args.n_comps + int(args.drop_first), min(matrix.shape) - 1)
    if n_comps < 1 + int(args.drop_first):
        raise SystemExit(f"cannot compute LSI of a {matrix.shape[0]} x {matrix.shape[1]} matrix")
    if n_comps < args.n_comps + int(args.drop_first):
        scio.log.warning("%d components requested, %d possible", args.n_comps, n_comps - int(args.drop_first))

    if args.gpu:
        scio.gpu_init()
        import cupy as xp
        from cupyx.scipy.sparse.linalg import svds

        u, s, vt = svds(scio.to_gpu(matrix), k=n_comps)
    else:
        import numpy as xp
        from scipy.sparse.linalg import svds

        u, s, vt = svds(matrix.astype(np.float64), k=n_comps)
    order = xp.argsort(s)[::-1]
    u, s, vt = u[:, order], s[order], vt[order]
    signs = xp.sign(u[xp.argmax(xp.abs(u), axis=0), xp.arange(u.shape[1])])
    signs[signs == 0] = 1
    u, vt = u * signs, vt * signs[:, None]
    u, s, vt = (scio.to_cpu(a) for a in (u, s, vt))
    if args.drop_first:
        u, s, vt = u[:, 1:], s[1:], vt[1:]
    if not args.no_scale_embeddings:
        u = (u - u.mean(axis=0)) / u.std(axis=0)

    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("obsm/X_lsi", u.astype(np.float32))
    writer.add("varm/LSI", vt.T.astype(np.float32))
    writer.add("uns/lsi", {"stdev": s / np.sqrt(reader.n_obs - 1)})
    writer.write(
        params={"matrix_key": args.matrix_key, "n_comps": int(u.shape[1]), "scale_embeddings": not args.no_scale_embeddings,
                "drop_first": args.drop_first, "gpu": args.gpu},
        versions=scio.versions("anndata", "numpy", "scipy", "cupy-cuda12x"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
