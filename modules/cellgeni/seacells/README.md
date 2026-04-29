# cellgeni/seacells

## Summary

Runs the SEACells algorithm to aggregate single-cell profiles into a smaller set of representative metacells (cellular states) from single-cell genomics data (e.g. scRNA-seq or scATAC-seq).

## Get started

Include this module in your Nextflow pipeline:

```nextflow
include { SEACELLS } from 'cellgeni/seacells'
```

### Inputs

- `tuple val(meta), path(adata)`
  - `meta`: sample metadata map (expects an `id` key; used for tagging and default output prefixing)
  - `adata`: input AnnData object in `.h5ad` format

### Outputs

- `h5ad`: `tuple val(meta), path("*.h5ad")` (metacell AnnData object; typically includes `seacell_metacells.h5ad`)
- `csv`: `tuple val(meta), path("*.csv")` (assignment tables; hard/soft assignments)
- `npy`: `tuple val(meta), path("*.npy")` (assignment weights)
- `pdf`: `tuple val(meta), path("*.pdf")` (diagnostic/QC plots)
- `pkl`: `tuple val(meta), path("*.pkl")` (pickled SEACells model)
- `versions`: `path("versions.yml")` (tool/package versions captured at runtime)

### Parameters

This module supports passing arguments to `seacells_aggregate.py` via `task.ext.args`.

When running the module directly with `nextflow module run`, you set these at the command line using Nextflow “process options”:

- `-process.ext.args='<SEACELLS_ARGS>'` (i.e. `-process.ext.args=<SEACELLS_ARGS>`)
- `-process.ext.prefix='<SAMPLE_PREFIX>'` (i.e. `-process.ext.prefix=<SAMPLE_PREFIX>`, optional; defaults to `${meta.id}`)

Defaults are:

`--type gex --n_top_genes 2000 --n_components 50 --convergence_epsilon 0.00001 --min_iter 10 --max_iter 50`

The output/sample prefix can be controlled via `task.ext.prefix` (defaults to `${meta.id}`).

#### SEACells arguments

These are the supported arguments you can include in `ext.args`:

- `--type` (required): input data modality. Choices: `gex` or `atac`.
- `--n_metacells` (optional): explicit number of metacells to infer.
- `--gamma` (optional): sets metacell count adaptively as $n_{metacells} = \mathrm{round}(n_{cells} / \gamma)$. Mutually exclusive with `--n_metacells`.
- `--n_top_genes` (optional, default `2000`): number of highly-variable genes to select (GEX preprocessing).
- `--n_components` (optional, default `50`): number of components used for PCA (GEX) or LSI (ATAC).
- `--convergence_epsilon` (optional, default `1e-5`): convergence threshold for optimization.
- `--min_iter` (optional, default `10`): minimum iterations for the SEACells fitting step.
- `--max_iter` (optional, default `50`): maximum iterations for the SEACells fitting step.
- `--n_waypoint_eigs` (optional, default `10`): number of waypoint eigenvectors used during initialization.
- `--use_sparse` (flag): use sparse matrix operations (if supported by the input representation).
- `--precomputed` (optional): key in `adata.obsm` containing a precomputed embedding to use instead of running PCA/LSI.
- `--celltype_label` (optional): obs column name used to compute cell type purity (QC plot).
- `--delimiter` (optional): append sample suffix to barcodes using this delimiter.

Notes:

- `--adata`, `--sample`, and `--output_dir` are handled by the module wrapper and do not need to be provided in `ext.args`.
- You must provide exactly one of `--n_metacells` or `--gamma`.

#### Full `nextflow module run` example

```bash
/software/cellgen/cellgeni/nextflow/26.04.0/nextflow module run cellgeni/seacells \
  --meta.id pbmc_10k \
  --adata /path/to/pbmc10k.h5ad \
  -process.ext.prefix=pbmc_10k \
  -process.ext.args='--type gex --gamma 75 --n_top_genes 2000 --n_components 50 --convergence_epsilon 0.00001 --min_iter 10 --max_iter 50'
```

## Dependencies

This module runs inside the container `quay.io/cellgeni/seacells:latest`.

## Citation

Persad S, Choo Z-N, Dien C, et al. *SEACells: Inference of transcriptional and epigenomic cellular states from single-cell genomics data*. bioRxiv (2022). https://doi.org/10.1101/2022.04.02.486748

## License

MIT
