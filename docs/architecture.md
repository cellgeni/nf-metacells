# nf-metacells architecture

Status: redesign in progress (started November 2026). Tag `0.1.0` is the last legacy release. This page
describes the target design that new code must follow; the roadmap at the end tracks what exists.

## Why the redesign

The legacy `SEACELLS` and `HMETACELLS` modules each ran their own preprocessing inside one Python script
(normalize → log1p → HVG → PCA, or TF-IDF → LSI), CPU only, and wrote whole AnnData copies. The redesign splits
every processing step into its own module with a CPU (scanpy/snapATAC2) and a GPU (rapids-singlecell, cuML,
CuPy) path, keeps intermediate outputs small, and separates *which cells form a metacell* (membership) from
*what gets summed* (aggregation), so new metacell methods are cheap to add.

## Layers

```
samplesheet: one row per sample (id, rna?, atac?, metadata?, optional precomputed step outputs)
 ├─ INGEST     adata-cli               ADATA_CONVERT → ADATA_ATTACHOBS                               CPU
 ├─ REPRESENT  only the steps that     RNA : NORMALIZETOTAL → LOG1P → HVG → PCA                      CPU|GPU
 │             the methods need        ATAC: TFIDF → LSI      and/or   SELECTFEATURES → SPECTRAL     CPU|GPU
 │                                     NEIGHBORS on the chosen embedding when a method needs a graph CPU|GPU
 ├─ CONSTRUCT  one task per method×run SEACELLS, HMETACELLS (more to come)                           CPU|GPU
 │                                     → each emits a membership table and nothing else is required
 ├─ AGGREGATE  method-agnostic         membership × every assay of the sample → metacell AnnData     CPU|GPU
 ├─ EVALUATE   method-agnostic         size, coverage, compactness, separation, purity               CPU
 └─ ASSEMBLE   on by default           base object + step outputs + membership columns (one copy)    CPU
```

Membership comes from one assay (`--membership_assay rna|atac`; joint later) and is applied to every assay of
the sample, so paired multiome cells are pooled consistently.

## Modules

Source lives in `modules/cellgeni/<tool>/<step>/` and is published to the `cellgeni` namespace on
registry.nextflow.io.

| Module | CPU | GPU (`task.accelerator`) | Output field |
|---|---|---|---|
| adata/convert | `adata convert` | – | `adata` |
| adata/attachobs | adata-cli export → duckdb join → subset → import | – | `adata` |
| scanpy/normalizetotal | `sc.pp.normalize_total` | `rsc.pp.normalize_total` | `normalized` |
| scanpy/log1p | `sc.pp.log1p` | `rsc.pp.log1p` | `log1p` |
| scanpy/hvg | `sc.pp.highly_variable_genes` (mask only) | `rsc.pp.highly_variable_genes` | `hvg` |
| scanpy/pca | `sc.pp.pca` on the HVG mask (`--scale` scales those genes in memory first) | `rsc.pp.pca` | `pca` |
| scanpy/neighbors | `sc.pp.neighbors` | `rsc.pp.neighbors` | `neighbors` |
| atac/tfidf | muon's formula, idf = 0 for empty peaks | same code on CuPy | `tfidf` |
| atac/lsi | sparse SVD + muon's ordering and scaling | CuPy sparse SVD | `lsi` |
| snapatac2/selectfeatures | `snap.pp.select_features` (mask only) | – | `selected` |
| snapatac2/spectral | `snap.tl.spectral` | CuPy port of the same algorithm | `spectral` |
| seacells | SEACells dense or sparse | SEACells unified GPU model | `membership` |
| hmetacells | per-group kNN + sknetwork hierarchy | per-group kNN on GPU | `membership` |
| scanpy/aggregate | `sc.get.aggregate` | `rsc.get.aggregate` | `metacells` |
| metacells/evaluate | QC metrics | – | `metrics` |
| anndata/assemble | copy base once, write elements | – | `final` |

## Records

```groovy
record Sample     { id; rna?; atac?; metadata?; embedding_key?; <one optional field per step> }  // samplesheet row
record AssayState { id; sample; assay; adata; counts_key; embedding_key?; todo; <one field per step> }
record Membership { id; sample; basis; method; run; membership; extras }
```

- Every channel item has a unique `id`: `S1.rna` for an assay state, `S1.rna:pca.seacells.r1`-style ids for
  memberships. Joins between different kinds of items use `sample`.
- A step filters the states that need it, runs, and joins its output back by `id`. `todo` (the steps a state
  still needs) is computed once from the method-needs table, what the samplesheet already provides and
  `--embedding_key`.
- Modules declare their record input inline; pipeline types live in `workflows/types.nf`.

## I/O contract

| `--output_mode` | What a step writes |
|---|---|
| `raw` (default) | `<id>.<step>.raw/` laid out like the AnnData store (`obsm/X_pca.npy`, `varm/PCs.npy`, `uns/pca.json`, `layers/log1p.npz`, `var/…csv`, `obsp/connectivities.npz`) plus `manifest.json` |
| `slim` | `<id>.<step>.h5ad\|.zarr` with obs/var names and only the new elements |
| `attach` | `<id>.<step>.h5ad\|.zarr` = input object + new elements (one full copy per step) |

- Formats: dense `.npy`; sparse `.npz` (uncompressed) or `.mtx` with `--sparse_format mtx`; tables `.csv`;
  dicts `.json`; objects `.h5ad` or `.zarr` (`--adata_format`).
- `manifest.json` records n_obs, n_vars and digests of the obs/var names; loaders refuse misaligned inputs.
- Any element input may be a raw directory, a slim object or a full object. Steps never change n_obs or n_vars
  (feature selection is a mask) and never overwrite the counts.
- The output index uses the samplesheet columns, so a finished run can be the next run's `--input`.
- Step scripts implement this contract through `shared/python/scio.py`, which `scripts/sync_shared.sh` copies into
  each module: `Reader(base, sources)` finds an element in the given step outputs, then in the base object, and
  refuses sources whose cells or features differ from the base; `Writer(reader, step, prefix, mode)` writes the new
  elements in any mode. Slim and attach objects keep the manifest in `uns['nf_metacells'][<step>]`.

## CPU and GPU

- Modules switch container and code path on `task.accelerator`; they never set it.
- `--gpu_steps none|auto|all|<list>` sets `accelerator` per process; LSF ignores `accelerator`, so the
  pipeline config adds `-gpu "mode=shared:j_exclusive=yes:gmem=…:num=1"`, the GPU queue and Singularity `--nv`.
- Farm constraints: driver 550 (CUDA 12.4) → CUDA 12 builds only; rapids-singlecell ≥ 0.15 lacks V100 kernels
  → `gpu_gmem` ≥ 40000 MB so jobs land on A100/H100.

## Adding a metacell method

1. A module `cellgeni/<method>` taking the inputs it needs (embedding, graph, counts, labels) and writing
   `<id>.membership.csv` (`barcode,metacell`; empty metacell = unassigned) plus optional extras.
2. A row in the method-needs table and a small adapter in `subworkflows/local/construct/`.
3. A `<method>_*` parameter group in `nextflow_schema.json` and `params {}`.
4. A test run on `tests/data` and a docs page. Aggregation, evaluation, assembly and publishing need no change.

Labels have two roles: `--constrain_by` (build metacells within groups) and `--eval_labels` (purity only).
`--metacell_size` sets a shared target granularity so methods can be compared fairly.

## Roadmap

| Phase | Scope | State |
|---|---|---|
| 0 | Legacy close-out (seacells 1.1.0, tag 0.1.0), Nextflow 26.10.7, test data, spikes, CI | done |
| 1 | adata-cli ingestion (`adata/convert`, `adata/attachobs`) | done |
| 2 | Typed pipeline skeleton, nf-schema, output block | done: the v1 SEACells/hmetacells modules still preprocess internally |
| 3 | Containers, GPU plumbing, shared I/O helper | done: `scanpy` and `snapatac2` images wait for quay.io push access |
| 4 | Representation modules (CPU + GPU) | in progress: modules and wiring done (the v1 methods read the computed embedding); GPU parity on the farm and real-data timings pending |
| 5 | SEACells/hmetacells 2.0, aggregate, evaluate, assemble | – |
| 6 | Registry publishing, docs, first calendar release | – |
| 7 | New methods: k-means, MetaQ, SuperCell, MC2, mcRigor, EpiCarousel, LSMetacell | – |
