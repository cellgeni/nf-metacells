#!/usr/bin/env python3
"""Mark informative peaks or tiles in var/selected with snapATAC2 select_features; never subsets.

Features are ranked by total accessibility after the lower and upper quantile filters; features without counts are
never selected. --whitelist/--blacklist take BED files and need region names (chr:start-end) as var names.
Counts come from --counts-key (X, layers/<name> or raw/X) of the base object.
"""

import argparse
import sys

import scio

STEP = "selectfeatures"


def main(argv: list[str] | None = None) -> int:
    parser = scio.add_io_arguments(argparse.ArgumentParser(description=__doc__))
    parser.add_argument("--counts-key", default="X", help="element with the counts: X, layers/<name> or raw/X")
    parser.add_argument("--n-features", type=int, default=500000)
    parser.add_argument("--filter-lower-quantile", type=float, default=0.005)
    parser.add_argument("--filter-upper-quantile", type=float, default=0.005)
    parser.add_argument("--whitelist", help="BED file of regions to keep")
    parser.add_argument("--blacklist", help="BED file of regions to drop")
    parser.add_argument("--max-iter", type=int, default=1, help="iterations of feature selection (1: by counts)")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    scio.setup_logging()
    if args.gpu:
        scio.log.warning("selectfeatures has no GPU path: running on CPU")

    import anndata as ad
    import pandas as pd
    import snapatac2 as snap

    reader = scio.Reader(args.adata)
    adata = ad.AnnData(X=reader.read(args.counts_key), var=pd.DataFrame(index=pd.Index(reader.var_names)))
    options = {
        "n_features": args.n_features,
        "filter_lower_quantile": args.filter_lower_quantile,
        "filter_upper_quantile": args.filter_upper_quantile,
        "whitelist": args.whitelist,
        "blacklist": args.blacklist,
        "max_iter": args.max_iter,
    }
    snap.pp.select_features(adata, inplace=True, n_jobs=args.threads, verbose=False, **options)
    selected = int(adata.var["selected"].sum())
    if selected < 2:
        raise SystemExit(f"only {selected} features selected: lower the quantile filters or raise --n-features")

    writer = scio.writer_from_args(reader, STEP, args)
    writer.add("var/selected", adata.var["selected"])
    writer.add("uns/selectfeatures", {"n_selected": selected})
    writer.write(params={"counts_key": args.counts_key, **options},
                 versions=scio.versions("anndata", "snapatac2"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
