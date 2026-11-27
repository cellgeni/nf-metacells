#!/usr/bin/env python3
"""Assemble a base AnnData object and step outputs into one object (one copy of the base).

Every element of each --sources step output (raw directory, slim or full object) is written into a copy of the base
object; elements of an object source that the base already has are left as they are in the base. Later sources win
over earlier ones. The step manifests are kept in uns['nf_metacells'].
"""

import argparse
import json
import sys

import scio

STEP = "assemble"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", required=True, help="base AnnData object (.h5ad or .zarr)")
    parser.add_argument("--sources", nargs="*", default=[], help="step outputs to add")
    parser.add_argument("--prefix", required=True, help="output name prefix, usually the record id")
    parser.add_argument("--adata-format", choices=scio.ADATA_FORMATS, default="h5ad")
    args = parser.parse_args(argv)
    scio.setup_logging()

    reader = scio.Reader(args.adata, args.sources)
    base_keys = set(reader.base.keys())
    elements: dict[str, tuple] = {}
    manifests: dict[str, str] = {}
    for source in reader.sources:  # later sources win
        if isinstance(source, scio.RawStore):
            keys = source.keys()
            manifests[source.manifest["step"]] = json.dumps(source.manifest)
        else:
            keys = [k for k in source.keys() if k not in base_keys and k not in ("X", "raw/X")]
            if f"uns/{scio.UNS_KEY}" in keys:
                keys.remove(f"uns/{scio.UNS_KEY}")
                manifests.update(source.read(f"uns/{scio.UNS_KEY}"))
        for key in keys:
            elements[key] = (source, key)

    writer = scio.Writer(reader, STEP, args.prefix, mode="attach", adata_format=args.adata_format)
    for key, (source, _) in elements.items():
        writer.add(key, source.read(key))
    if manifests:  # the attach writer adds this step's manifest to them
        writer.add(f"uns/{scio.UNS_KEY}", manifests)
    writer.write(params={"sources": [str(s.path) for s in reader.sources]}, versions=scio.versions("anndata"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
