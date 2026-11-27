# cellgeni/atac/tfidf

TF-IDF of peak counts into layers/tfidf, with muon's formula (the v1 modules' ATAC preprocessing): NumPy/SciPy on CPU, CuPy when the task has an accelerator.

```nextflow
nextflow.enable.types = true
include { ATAC_TFIDF } from 'cellgeni/atac/tfidf'

ATAC_TFIDF(channel.of(record(id: 'S1.atac', adata: file('S1_atac.h5ad'), counts_key: null, output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `counts_key` | String? | Optional element with the counts: `X` (default), `layers/<name>` or `raw/X`. |
| `output_prefix` | String? | Output prefix (default `id`) |

Step outputs given as inputs (raw directories, slim or full objects) must hold the same cells and features as the base
object, in the same order; elements missing from them are read from the base object.

## Output record

`record(id, tfidf)`, where `tfidf` is `<prefix>.tfidf.raw` or `<prefix>.tfidf.h5ad|zarr`, holding
layers/tfidf.

## Options

- `ext.args`: `--scale-factor <float>` (default 1e4), `--no-log-tf`, `--no-log-idf`, `--log-tfidf` (with both previous). Peaks and cells without counts stay zero; no peak is removed.
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90` (CPU: scanpy 1.12.6,
anndata 0.13.4), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are emitted on the `versions` topic.
