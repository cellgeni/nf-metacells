# cellgeni/anndata/assemble

Assembles a base AnnData object and step outputs into one object: a copy of the base with every element of the step
outputs written into it. Step outputs can be raw step directories (`<prefix>.<step>.raw`, see
`cellgeni/scanpy/*`), slim objects or full objects; object sources contribute only the elements the base lacks.
Later sources win. The step manifests are kept in `uns['nf_metacells']`.

```nextflow
nextflow.enable.types = true
include { ANNDATA_ASSEMBLE } from 'cellgeni/anndata/assemble'

ANNDATA_ASSEMBLE(channel.of(record(id: 'S1.rna', adata: file('S1.h5ad'), sources: [file('S1.rna.pca.raw')], output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/` |
| `sources` | List\<Path\> | Step outputs to add, in order; they must hold the same cells and features as the base |
| `output_prefix` | String? | Output prefix (default `id`) |

## Output record

`record(id, assembled)`, where `assembled` is `<prefix>.assemble.h5ad` (or `.zarr`).

## Options

- `ext.adata_format`: `h5ad` (default) or `zarr`.

Container: `quay.io/cellgeni/scanpy:1.12.6`. Versions are emitted on the `versions` topic.
