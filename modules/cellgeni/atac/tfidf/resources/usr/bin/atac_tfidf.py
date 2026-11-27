#!/usr/bin/env python3
"""TF-IDF of peak counts into layers/tfidf, with muon's formula (the v1 modules' ATAC preprocessing).

TF = counts / total counts of the cell x --scale-factor, then log1p (unless --no-log-tf); IDF = number of cells /
total counts of the peak, then log1p (unless --no-log-idf); TF-IDF = TF x IDF, then log1p with --log-tfidf (which
needs --no-log-tf and --no-log-idf). Peaks and cells without counts stay zero; no peak is removed. NumPy/SciPy on
CPU, CuPy with --gpu. Counts come from --counts-key (X, layers/<name> or raw/X) of the base object.
"""

import argparse
import sys

import scio

STEP = "tfidf"


def tfidf(counts, scale_factor, log_tf, log_idf, log_tfidf, gpu=False):
    """CSR float32 TF-IDF of a CSR count matrix (scipy, or cupyx with gpu=True)."""
    if gpu:
        import cupy as xp
        import cupyx.scipy.sparse as sparse
    else:
        import numpy as xp
        import scipy.sparse as sparse
    counts = counts.tocsr().astype(xp.float32)
    n_cells = counts.shape[0]
    cell_totals = xp.asarray(counts.sum(axis=1)).ravel()
    peak_totals = xp.asarray(counts.sum(axis=0)).ravel()
    # Empty cells and peaks: zero instead of muon's inf x 0 = nan -> 0
    inv_cell = xp.where(cell_totals > 0, 1.0 / xp.maximum(cell_totals, 1), 0).astype(xp.float32)
    tf = (sparse.diags(inv_cell) @ counts).tocsr()
    if scale_factor not in (None, 0, 1):
        tf.data *= xp.float32(scale_factor)
    if log_tf:
        tf.data = xp.log1p(tf.data)
    idf = xp.where(peak_totals > 0, n_cells / xp.maximum(peak_totals, 1), 0).astype(xp.float32)
    if log_idf:
        idf = xp.log1p(idf)
    result = (tf @ sparse.diags(idf)).tocsr()
    if log_tfidf:
        result.data = xp.log1p(result.data)
    result.eliminate_zeros()
    return result


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--counts-key", default="X", help="element with the peak counts: X, layers/<name> or raw/X")
    parser.add_argument("--scale-factor", type=float, default=1e4, help="multiplies TF before the log (default 1e4)")
    parser.add_argument("--no-log-tf", action="store_true")
    parser.add_argument("--no-log-idf", action="store_true")
    parser.add_argument("--log-tfidf", action="store_true", help="log1p of TF x IDF (with --no-log-tf --no-log-idf)")
    args = parser.parse_args(argv)
    scio.setup_logging()
    if args.log_tfidf and not (args.no_log_tf and args.no_log_idf):
        raise SystemExit("--log-tfidf needs --no-log-tf and --no-log-idf (as in muon)")

    reader = scio.Reader(args.adata)
    counts = reader.read(args.counts_key)
    options = {"scale_factor": args.scale_factor, "log_tf": not args.no_log_tf, "log_idf": not args.no_log_idf,
               "log_tfidf": args.log_tfidf}
    if args.gpu:
        scio.gpu_init()
        counts = scio.to_gpu(counts)
    result = tfidf(counts, gpu=args.gpu, **options)

    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("layers/tfidf", result)
    writer.write(params={"counts_key": args.counts_key, "gpu": args.gpu, **options},
                 versions=scio.versions("anndata", "numpy", "scipy", "cupy-cuda12x"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
