# cellgeni/scanpy/neighbors

k-nearest-neighbour graph of an embedding into obsp/connectivities, obsp/distances and uns/neighbors: scanpy on CPU, rapids-singlecell when the task has an accelerator.

```nextflow
nextflow.enable.types = true
include { SCANPY_NEIGHBORS } from 'cellgeni/scanpy/neighbors'

SCANPY_NEIGHBORS(channel.of(record(id: 'S1.rna', adata: file('S1.h5ad'), embedding: null, embedding_key: null, output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `embedding` | Path? | Optional step output with the embedding (default: the base object). |
| `embedding_key` | String? | Optional embedding element (default `obsm/X_pca`). |
| `output_prefix` | String? | Output prefix (default `id`) |

Step outputs given as inputs (raw directories, slim or full objects) must hold the same cells and features as the base
object, in the same order; elements missing from them are read from the base object.

## Output record

`record(id, neighbors)`, where `neighbors` is `<prefix>.neighbors.raw` or `<prefix>.neighbors.h5ad|zarr`, holding
obsp/connectivities, obsp/distances, uns/neighbors.

## Options

- `ext.args`: `--n-neighbors <int>` (default 15), `--n-dims <int>` (first columns only), `--metric` (default euclidean), `--method umap|gauss`, `--seed` (default 0). The GPU search is exact; scanpy's is approximate on large data.
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90` (CPU: scanpy 1.12.6,
anndata 0.13.4), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are emitted on the `versions` topic.
