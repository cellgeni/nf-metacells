# cellgeni/atac/lsi

Latent semantic indexing of a TF-IDF layer into obsm/X_lsi, varm/LSI and uns/lsi, as muon's lsi does: SciPy svds on CPU, cupyx svds when the task has an accelerator.

```nextflow
nextflow.enable.types = true
include { ATAC_LSI } from 'cellgeni/atac/lsi'

ATAC_LSI(channel.of(record(id: 'S1.atac', adata: file('S1_atac.h5ad'), tfidf: null, output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `tfidf` | Path? | Optional step output with the TF-IDF layer (default: the base object). |
| `output_prefix` | String? | Output prefix (default `id`) |

Step outputs given as inputs (raw directories, slim or full objects) must hold the same cells and features as the base
object, in the same order; elements missing from them are read from the base object.

## Output record

`record(id, lsi)`, where `lsi` is `<prefix>.lsi.raw` or `<prefix>.lsi.h5ad|zarr`, holding
obsm/X_lsi, varm/LSI, uns/lsi.

## Options

- `ext.args`: `--n-comps <int>` (default 50), `--matrix-key` (default `layers/tfidf`), `--no-scale-embeddings`, `--drop-first`. Components are signed with sklearn's svd_flip convention.
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions). The
  GPU path needs an A100 or H100 (CUDA 12).

Containers: `community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90` (CPU: scanpy 1.12.6,
anndata 0.13.4), `quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are emitted on the `versions` topic.
