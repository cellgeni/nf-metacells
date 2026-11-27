# cellgeni/scanpy/log1p

Log-transforms a normalised layer, log(1 + x), into layers/log1p: scanpy on CPU, rapids-singlecell when the task has an accelerator.

```nextflow
nextflow.enable.types = true
include { SCANPY_LOG1P } from 'cellgeni/scanpy/log1p'

SCANPY_LOG1P(channel.of(record(id: 'S1.rna', adata: file('S1.h5ad'), normalized: null, output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `normalized` | Path? | Optional step output with the normalised layer (default: the base object). |
| `output_prefix` | String? | Output prefix (default `id`) |

Step outputs given as inputs (raw directories, slim or full objects) must hold the same cells and features as the base
object, in the same order; elements missing from them are read from the base object.

## Output record

`record(id, log1p)`, where `log1p` is `<prefix>.log1p.raw` or `<prefix>.log1p.h5ad|zarr`, holding
layers/log1p.

## Options

- `ext.args`: `--layer-key <element>` (default `layers/normalized`), `--base <float>` (default: natural logarithm).
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `quay.io/cellgeni/scanpy:1.12.6` (CPU), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are
emitted on the `versions` topic.
