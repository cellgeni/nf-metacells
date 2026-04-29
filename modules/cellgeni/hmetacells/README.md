# cellgeni/hmetacells

## Summary

Runs the hierarchical metacell aggregation workflow to cluster single-cell profiles into a smaller set of representative metacells from single-cell genomics data.

## Get started

Include this module in your Nextflow pipeline:

```nextflow
include { HMETACELLS } from 'cellgeni/hmetacells'
```

### Inputs

- `tuple val(meta), path(adata)`
  - `meta`: sample metadata map (expects an `id` key; used for tagging and default output prefixing)
  - `adata`: input AnnData object in `.h5ad` format

### Outputs

- `csv`: `tuple val(meta), path("*.csv")` (hierarchical metacell assignment table; typically includes `hierarchial_metacells.csv`)
- `versions`: `path("versions.yml")` (tool/package versions captured at runtime)

### Parameters

This module supports passing arguments to `hierarchial_metacells.py` via `task.ext.args`.

When running the module directly with `nextflow module run`, set these at the command line using Nextflow “process options”:

- `-process.ext.args='<HMETACELLS_ARGS>'` (i.e. `-process.ext.args=<HMETACELLS_ARGS>`)
- `-process.ext.prefix='<SAMPLE_PREFIX>'` (i.e. `-process.ext.prefix=<SAMPLE_PREFIX>`, optional; defaults to `${meta.id}`)

Defaults are:

`--type gex --n_min 5 --n_max 20 --method ward --n_top_genes 2000 --n_components 50 --n_neighbors 15`

The output/sample prefix can be controlled via `task.ext.prefix` (defaults to `${meta.id}`).

#### Hierarchical aggregation arguments

These are the supported arguments you can include in `ext.args`:

- `--type` (required): input data modality. Choices: `gex` or `atac`.
- `--n_min` (required): minimum number of cells in a metacell.
- `--n_max` (required): maximum number of cells in a metacell.
- `--method` (required): clustering method. Choices: `paris`, `louvain`, or the default workflow value.
- `--n_top_genes` (optional, default `2000`): number of highly-variable genes to select for GEX preprocessing.
- `--n_components` (optional, default `50`): number of components used for PCA (GEX) or LSI (ATAC).
- `--n_neighbors` (optional, default `15`): number of nearest neighbors used during clustering.
- `--celltype_label` (optional): obs column used to split cells by cell type before clustering.
- `--precomputed` (optional): key in `adata.obsm` containing a precomputed embedding.
- `--delimiter` (optional): append sample suffix to barcodes using this delimiter.

Notes:

- `--adata`, `--sample`, and `--output` are handled by the module wrapper and do not need to be provided in `ext.args`.

#### Full `nextflow module run` example

```bash
/software/cellgen/cellgeni/nextflow/26.04.0/nextflow module run cellgeni/hmetacells \
  --meta.id pbmc_10k \
  --adata /path/to/pbmc10k.h5ad \
  -process.ext.prefix=pbmc_10k \
  -process.ext.args='--type gex --n_min 5 --n_max 20 --method louvain --n_top_genes 2000 --n_components 50 --n_neighbors 15 --celltype_label celltype'
```

## Dependencies

This module runs inside the container `quay.io/cellgeni/metacells-python:latest`.

## Citation

Add the relevant citation for the hierarchical metacell method here.

## License

MIT
