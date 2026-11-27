# cellgeni/snapatac2/spectral

Spectral embedding (Laplacian eigenmaps of cosine similarity) of peak counts into obsm/X_spectral and uns/spectral_eigenvalue: snapATAC2 on CPU, a CuPy port of its matrix-free cosine algorithm when the task has an accelerator.

```nextflow
nextflow.enable.types = true
include { SNAPATAC2_SPECTRAL } from 'cellgeni/snapatac2/spectral'

SNAPATAC2_SPECTRAL(channel.of(record(id: 'S1.atac', adata: file('S1_atac.h5ad'), selected: null, counts_key: null, output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `selected` | Path? | Optional step output with the feature mask (default: the base object). |
| `counts_key` | String? | Optional element with the counts: `X` (default), `layers/<name>` or `raw/X`. |
| `output_prefix` | String? | Output prefix (default `id`) |

Step outputs given as inputs (raw directories, slim or full objects) must hold the same cells and features as the base
object, in the same order; elements missing from them are read from the base object.

## Output record

`record(id, spectral)`, where `spectral` is `<prefix>.spectral.raw` or `<prefix>.spectral.h5ad|zarr`, holding
obsm/X_spectral, uns/spectral_eigenvalue.

## Options

- `ext.args`: `--n-comps <int>` (default 30; fewer columns when non-positive eigenvalues are dropped), `--features-key` (default `var/selected`), `--all-features`, `--distance-metric cosine|jaccard`, `--sample-size` (Nystrom, CPU), `--no-weighted-by-sd`, `--seed` (default 0). The GPU path computes the full cosine embedding only.
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `quay.io/cellgeni/snapatac2:2.10.0` (CPU), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are emitted on the `versions` topic.
