# cellgeni/hmetacells

Hierarchical metacells built within the groups of an obs column. Each group's kNN graph of a precomputed embedding
is cut by a scikit-network hierarchy into balanced clusters of at most `--n-max` cells.

```nextflow
nextflow.enable.types = true
include { HMETACELLS } from 'cellgeni/hmetacells'

HMETACELLS(channel.of(record(id: 'S1.rna', adata: file('S1_rna.h5ad'), embedding: file('S1.pca.raw'),
                             embedding_key: 'obsm/X_pca', output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/`; supplies the barcodes and the grouping column |
| `embedding` | Path? | Step output holding the embedding (default: the base object) |
| `embedding_key` | String? | Element key of the embedding (default `obsm/X_pca`) |
| `output_prefix` | String? | Output prefix (default `id`) |

## Output record

`record(id, membership, summary)`:

- `<prefix>.membership.csv` — `barcode,metacell`, one row per cell. Metacells are named `<group>_<n>`; a cell
  whose group is missing gets an empty `metacell`, which means unassigned.
- `<prefix>.hmetacells.json` — `n_cells`, `n_metacells`, `n_groups`, `unassigned`, `metacell_size`
  (min/median/max), `below_n_min`, `constrain_by`, `embedding_key`, `hierarchy`, `n_neighbors`, `n_min`,
  `n_max`, `gpu`, `seed` and tool versions.

## Options

- `ext.args`: `--constrain-by <obs column>`, `--n-min <int>` and `--n-max <int>` are **required**. Also
  `--method louvain|paris` (default `louvain`), `--n-neighbors` (default 15), `--n-dims` (use the first n columns
  of the embedding) and `--seed` (default 0).
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions).
  Only the kNN graphs move to the device; the hierarchy and the cut stay on the CPU. The GPU path needs an A100
  or H100 (CUDA 12).

## Behaviour worth knowing

- **`--n-max` is a bound, `--n-min` is a report.** `cut_balanced` bounds metacells from above only, so a cluster
  below `--n-min` says that part of the group's graph is sparse. Those are counted in `below_n_min` and left
  alone rather than merged, which would undo the balance the cut just achieved.
- **A group smaller than `--n-max` becomes one metacell**, without building a graph.
- **The neighbourhood shrinks to fit a group.** `--n-neighbors` is capped at one less than the group's size,
  because both backends refuse to ask for more neighbours than there are other cells.

Containers: `community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90` (CPU),
`quay.io/cellgeni/rapids-singlecell:0.18.1` (GPU). Versions are emitted on the `versions` topic.
