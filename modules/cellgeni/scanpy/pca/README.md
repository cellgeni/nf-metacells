# cellgeni/scanpy/pca

Principal components of a log-normalised layer, on the highly variable genes, into obsm/X_pca, varm/PCs and uns/pca: scanpy on CPU, rapids-singlecell when the task has an accelerator.

```nextflow
nextflow.enable.types = true
include { SCANPY_PCA } from 'cellgeni/scanpy/pca'

SCANPY_PCA(channel.of(record(id: 'S1.rna', adata: file('S1.h5ad'), matrix: null, mask: null, output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `matrix` | Path? | Optional step output with the log-normalised layer (default: the base object). |
| `mask` | Path? | Optional step output with the gene mask (default: the base object). |
| `output_prefix` | String? | Output prefix (default `id`) |

Step outputs given as inputs (raw directories, slim or full objects) must hold the same cells and features as the base
object, in the same order; elements missing from them are read from the base object.

## Output record

`record(id, pca)`, where `pca` is `<prefix>.pca.raw` or `<prefix>.pca.h5ad|zarr`, holding
obsm/X_pca, varm/PCs (zero for unmasked genes), uns/pca.

## Options

- `ext.args`: `--n-comps <int>` (default 50; capped by the data), `--matrix-key` (default `layers/log1p`), `--mask-key` (default `var/highly_variable`), `--no-mask`, `--no-zero-center`, `--svd-solver`, `--scale` and `--scale-max-value` (scale the masked genes in memory first), `--seed` (default 0).
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90` (CPU: scanpy 1.12.6,
anndata 0.13.4), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are emitted on the `versions` topic.
