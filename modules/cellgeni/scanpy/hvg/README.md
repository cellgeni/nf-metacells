# cellgeni/scanpy/hvg

Flags highly variable genes in var/highly_variable, with the flavour's statistics in other var columns; never subsets: scanpy on CPU, rapids-singlecell when the task has an accelerator.

```nextflow
nextflow.enable.types = true
include { SCANPY_HVG } from 'cellgeni/scanpy/hvg'

SCANPY_HVG(channel.of(record(id: 'S1.rna', adata: file('S1.h5ad'), log1p: null, counts_key: null, output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `log1p` | Path? | Optional step output with the log-normalised layer (default: the base object). |
| `counts_key` | String? | Optional counts element for the seurat_v3 flavours: `X` (default), `layers/<name>` or `raw/X`. |
| `output_prefix` | String? | Output prefix (default `id`) |

Step outputs given as inputs (raw directories, slim or full objects) must hold the same cells and features as the base
object, in the same order; elements missing from them are read from the base object.

## Output record

`record(id, hvg)`, where `hvg` is `<prefix>.hvg.raw` or `<prefix>.hvg.h5ad|zarr`, holding
var/highly_variable and the flavour's statistics (var/means, var/dispersions_norm, ...), uns/hvg.

## Options

- `ext.args`: `--flavor seurat|cell_ranger|seurat_v3|seurat_v3_paper` (default seurat), `--n-top-genes <int>`, `--batch-key <obs column>`, `--min-mean`, `--max-mean`, `--min-disp`, `--span`, `--n-bins`, `--layer-key <element>` (default `layers/log1p`).
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `quay.io/cellgeni/scanpy:1.12.6` (CPU), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are
emitted on the `versions` topic.
