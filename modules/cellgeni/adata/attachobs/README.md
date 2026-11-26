# cellgeni/adata/attachobs

Attaches cell metadata to an AnnData store's `obs` and/or makes barcodes unique per sample, using
[adata-cli](https://github.com/cellgeni/adata-cli) and duckdb. Works on `.h5ad` and `.zarr` stores.

```nextflow
nextflow.enable.types = true
include { ADATA_ATTACHOBS } from 'cellgeni/adata/attachobs'

ADATA_ATTACHOBS(channel.of(record(
    id: 'S1', adata: file('S1.h5ad'), metadata: file('metadata.csv'),
    sample_id: 'S1', output_prefix: null, meta: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output name and sample id |
| `adata` | Path | `.h5ad` or `.zarr` store (staged under `input/`) |
| `metadata` | Path? | CSV with a `barcode` column plus annotations; a `sample_id` column restricts it to this sample |
| `sample_id` | String? | Selects the metadata rows when the CSV has a `sample_id` column; used for the barcode suffix (default `id`) |
| `output_prefix` | String? | Output name without extension (default `id`) |
| `meta` | Record? | Passed through unchanged |

## Output record

`record(id, meta, adata)` where `adata` is `<output_prefix or id>.h5ad` (or `.zarr`).

## Options (`ext.args`)

| Option | Effect |
|---|---|
| `--barcode-column COL` | obs column joined to the metadata `barcode` (default: obs names) |
| `--join inner\|left` | `inner` (default) drops cells without metadata; `left` keeps them with empty values |
| `--delimiter STR` | Rename obs names to `<barcode><STR><sample_id>`; originals go to a `barcode` column |
| `--sample-column NAME` | Add (or overwrite) a column holding the sample id |
| `--allow-dtype-change` | Accept categorical columns losing their type (see below) |

`ext.extension` sets the output format: `h5ad` (default) or `zarr`.

Behaviour: rows keep their original order; duplicated barcodes in the metadata are an error; metadata columns that
already exist in `obs` are attached with a `_metadata` suffix; the join always uses the original barcodes.

Known limitation: `obs` is round-tripped through CSV, and adata-cli 0.17.6 turns categorical columns with
digit-only labels (e.g. leiden `"0"`, `"1"`) into integers
([cellgeni/adata-cli#124](https://github.com/cellgeni/adata-cli/issues/124)). The module fails in that case unless
`--allow-dtype-change` is given.

Container: `quay.io/cellgeni/adata-cli:0.17.6` (adata-cli and duckdb). Versions are emitted on the `versions` topic.
