# nf-metacells

Nextflow pipeline that builds metacells from single-cell RNA and ATAC AnnData: [SEACells](https://github.com/dpeerlab/SEACells)
and hierarchical metacells built within groups (e.g. cell types). Inputs can be `.h5ad`/`.zarr` stores or raw 10x
output, and samples can have RNA, ATAC or both.

The pipeline is being redesigned into per-step CPU/GPU modules; see [docs/architecture.md](docs/architecture.md) for
the design and roadmap, and [CHANGELOG.md](CHANGELOG.md) for what has changed.

## Quick start (Sanger farm)

```bash
module load cellgen/nextflow/26.10.7
nextflow run cellgeni/nf-metacells -r dev --input samplesheet.csv --methods seacells --outdir results
```

The default configuration submits to LSF and runs the tools in Singularity containers. Add `-profile local` to run
everything on the current machine, or `-profile docker` to use Docker. `--help` lists the main parameters,
`--help <param>` explains one, `--helpFull` lists all of them.

## Samplesheet

One row per sample:

```csv
id,rna,atac,metadata
S1,/path/S1_gex.h5ad,/path/S1_atac.h5ad,
S2,/path/S2/filtered_feature_bc_matrix.h5,/path/S2/filtered_feature_bc_matrix.h5,/path/S2_metadata.csv
S3,/path/S3/filtered_feature_bc_matrix,,
```

- `id`: letters, digits, `_` or `-`.
- `rna`, `atac`: an `.h5ad` or `.zarr` store, a 10x `.h5` file or a 10x mtx folder; at least one per sample. 10x
  inputs are converted with [adata-cli](https://github.com/cellgeni/adata-cli). A multiome `.h5` can be given for both
  assays: each keeps its own feature type. A bare `.mtx` file is not accepted; give its folder.
- `metadata` (optional): a CSV with a `barcode` column and annotations, joined to that sample's cells.

Paths may be absolute or relative to the launch directory.

## Main parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | – | Samplesheet (required) |
| `--outdir` | `results` | Output directory |
| `--methods` | `seacells` | Comma-separated: `seacells`, `hmetacells` |
| `--membership_assay` | `rna` | Assay that drives membership in samples with both RNA and ATAC |
| `--metacell_size` | `75` | Target cells per metacell (SEACells); `--seacells_n_metacells` sets the number instead |
| `--constrain_by` | – | obs column to build metacells within (required by `hmetacells`) |
| `--eval_labels` | – | obs column used to report metacell purity (SEACells) |
| `--embedding_key` | – | obsm key of an existing embedding to use instead of computing one |
| `--atac_embedding` | `lsi` | ATAC embedding: `lsi` (TF-IDF + LSI, as muon) or `spectral` (snapATAC2) |
| `--hvg_n_top_genes`, `--pca_n_comps`, `--lsi_n_comps`, `--spectral_n_comps` | 2000, 50, 50, 30 | Embedding sizes |
| `--rna_counts_key`, `--atac_counts_key` | `X` | Where the counts are: `X`, `layers/<name>` or `raw/X` |
| `--gpu_steps` | `none` | Steps that run on a GPU: `none`, `auto`, `all` or a list (e.g. `pca,lsi,spectral`); A100/H100 only |
| `--cell_metadata` | – | CSV for all samples (`barcode`, `sample_id`, annotations); cells without metadata are dropped |
| `--barcode_column` | `obs_names` | obs column matched to the metadata barcodes |
| `--barcode_delimiter` | – | Rename barcodes to `<barcode><delimiter><id>` (originals kept in `obs['barcode']`, the id in `obs['sample']`) |
| `--allow_dtype_change` | false | Accept digit-only categorical obs columns (e.g. leiden) becoming integers when attaching metadata ([adata-cli#124](https://github.com/cellgeni/adata-cli/issues/124)); without it the run stops |

SEACells options start with `seacells_` (`n_metacells`, `n_waypoint_eigs`, `convergence_epsilon`, `min_iter`,
`max_iter`, `use_sparse`) and hierarchical-metacell options with `hmetacells_` (`n_min`, `n_max`, `method`,
`n_neighbors`); `--helpFull` describes them all.

Representation options (`normalize_target_sum`, `hvg_flavor`, `hvg_batch_key`, `pca_scale`, `tfidf_scale_factor`,
`lsi_drop_first`, `select_*`, `spectral_*`, `seed`) and the format of step outputs (`--output_mode raw|slim|attach`,
`--adata_format`, `--sparse_format`) are also in `--helpFull`.

Each method runs on one assay per sample (`--membership_assay`). The pipeline computes that assay's embedding one
step per task (RNA: normalise, log1p, highly variable genes, PCA; ATAC: TF-IDF and LSI, or feature selection and a
spectral embedding), each with a CPU and a GPU path; until the methods are rewritten (roadmap phase 5) they read the
embedding from a copy of the input object.

## Outputs

```
results/
├── <sample>/
│   ├── ingest/<sample>.<assay>.h5ad    AnnData the pipeline wrote (converted from 10x and/or with metadata attached)
│   ├── <assay>/<sample>.<assay>.<step>.raw   representation steps (normalizetotal, log1p, hvg, pca, tfidf, lsi,
│   │                                   selectfeatures, spectral): one file per element and manifest.json
│   ├── seacells/                       membership (seacell_metacell_assignments.csv), soft assignments,
│   │                                   metacell AnnData, model, plots
│   └── hmetacells/                     membership (hierarchial_metacells.csv)
├── memberships.csv                     one row per sample and method, with the membership table's path
├── ingested.csv                        one row per AnnData written under ingest/
├── representations.csv                 one row per representation step output
└── pipeline_info/                      run reports and software versions
```

## Test profiles

Synthetic data in [tests/data](tests/data); run from the repository root:

```bash
nextflow run . -profile test,local          # two RNA samples, SEACells and hmetacells
nextflow run . -profile test_atac,local     # one ATAC sample
nextflow run . -profile test_paired,local   # paired RNA+ATAC with cell metadata and barcode suffixes
nextflow run . -profile test_10x,local      # 10x multiome .h5 and mtx folder inputs
```

## Migrating from 0.1.0

| 0.1.0 | Now |
|---|---|
| `--filelist` (`id,path`) with `--type gex\|atac` | `--input` samplesheet with `rna` and/or `atac` columns |
| `--raw` | not needed: 10x inputs are detected and converted |
| `--seacells.enabled`, `--hierarchial.enabled` | `--methods seacells,hmetacells` |
| `--seacells.gamma` | `--metacell_size` |
| `--seacells.n_cells` | `--seacells_n_metacells` |
| `--seacells.n_top_genes`, `--hierarchial.n_top_genes` | `--hvg_n_top_genes` |
| `--seacells.n_components`, `--hierarchial.n_components` | `--pca_n_comps` (RNA), `--lsi_n_comps` (ATAC) |
| `--seacells.precomputed`, `--hierarchial.precomputed` | `--embedding_key` |
| `--seacells.min_iterations`, `--seacells.max_iterations` | `--seacells_min_iter`, `--seacells_max_iter` |
| `--seacells.<option>`, `--hierarchial.<option>` | `--seacells_<option>`, `--hmetacells_<option>` |
| `--celltype_label` | `--constrain_by` (hmetacells groups), `--eval_labels` (SEACells purity) |
| `--delimiter` | `--barcode_delimiter` |
| `results/<id>/…` | `results/<sample>/<method>/…` |

Nextflow stops with "Parameter … is not declared" when given a 0.1.0 parameter.
