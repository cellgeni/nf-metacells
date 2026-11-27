# cellgeni/snapatac2/selectfeatures

Marks informative peaks or tiles in var/selected with snapATAC2 select_features; never subsets: snapATAC2 on CPU.

```nextflow
nextflow.enable.types = true
include { SNAPATAC2_SELECTFEATURES } from 'cellgeni/snapatac2/selectfeatures'

SNAPATAC2_SELECTFEATURES(channel.of(record(id: 'S1.atac', adata: file('S1_atac.h5ad'), counts_key: null, output_prefix: null)))
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

`record(id, selected)`, where `selected` is `<prefix>.selectfeatures.raw` or `<prefix>.selectfeatures.h5ad|zarr`, holding
var/selected, uns/selectfeatures.

## Options

- `ext.args`: `--n-features <int>` (default 500000), `--filter-lower-quantile`, `--filter-upper-quantile` (default 0.005), `--whitelist <bed>`, `--blacklist <bed>` (need region var names), `--max-iter <int>` (default 1).
- `ext.output_mode`: `raw` (default; a directory laid out like the AnnData store, with `manifest.json`), `slim` (an
  object with only the new elements) or `attach` (a copy of the base object plus the new elements).
- `ext.adata_format`: `h5ad` (default) or `zarr`, for slim and attach outputs.
- `ext.sparse_format`: `npz` (default) or `mtx`, for raw outputs.
- CPU only; `--threads` follows `task.cpus`.

Containers: `quay.io/cellgeni/snapatac2:2.10.0`. Versions are emitted on the `versions` topic.
