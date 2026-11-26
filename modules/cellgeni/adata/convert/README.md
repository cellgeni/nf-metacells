# cellgeni/adata/convert

Converts a single-cell matrix to an AnnData store with [adata-cli](https://github.com/cellgeni/adata-cli):
10x `.h5` (including multiome), 10x mtx folders, h5ad, zarr, loom, MuData and more. Data are streamed from disk.

```nextflow
nextflow.enable.types = true
include { ADATA_CONVERT } from 'cellgeni/adata/convert'

ADATA_CONVERT(channel.of(record(
    id: 'S1', source: file('filtered_feature_bc_matrix.h5'),
    convert_args: '--feature-type "Gene Expression" --var-names gene_symbols --make-unique',
    output_prefix: null, meta: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output name |
| `source` | Path | Matrix to convert (staged under `input/`) |
| `convert_args` | String? | Per-item `adata convert` options, e.g. `--feature-type Peaks` for the ATAC half of a multiome `.h5` |
| `output_prefix` | String? | Output name without extension (default `id`) |
| `meta` | Record? | Passed through unchanged |

## Output record

`record(id, meta, adata)` where `adata` is `<output_prefix or id>.h5ad` (or `.zarr`).

## Options

- `ext.extension`: `h5ad` (default) or `zarr`.
- `ext.args`: extra `adata convert` options for every item.

Notes: features are named by gene id unless `--var-names gene_symbols` is given (add `--make-unique`); counts stay
integers; `--feature-type` only applies to 10x inputs; a bare `.mtx` file is refused — pass its 10x folder.

Container: `quay.io/cellgeni/adata-cli:0.17.6`. Versions are emitted on the `versions` topic.
