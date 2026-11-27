# cellgeni/scanpy/normalizetotal

Normalises counts to the same total per cell into `layers/normalized`: scanpy `normalize_total` on CPU,
rapids-singlecell when the task has an accelerator. Counts are cast to float32; the counts element is never changed.

```nextflow
nextflow.enable.types = true
include { SCANPY_NORMALIZETOTAL } from 'cellgeni/scanpy/normalizetotal'

SCANPY_NORMALIZETOTAL(channel.of(record(id: 'S1.rna', adata: file('S1.h5ad'), counts_key: 'layers/counts', output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `counts_key` | String? | Element with the counts: `X` (default), `layers/<name>` or `raw/X` |
| `output_prefix` | String? | Output prefix (default `id`) |

## Output record

`record(id, normalized)`, where `normalized` is `<prefix>.normalizetotal.raw` or `<prefix>.normalizetotal.h5ad|zarr`.

## Options

- `ext.args`: `--target-sum <float>` (default: the median total), `--exclude-highly-expressed`, `--max-fraction <float>`.
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the input plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90` (CPU: scanpy 1.12.6,
anndata 0.13.4), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are emitted on the `versions` topic.
