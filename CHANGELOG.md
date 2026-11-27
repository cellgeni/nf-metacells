# Changelog

All notable changes to nf-metacells. Pipeline releases use calendar versions (`YY-DDD`) from the redesign
onwards; modules under `modules/cellgeni/` use semver and have their own versions in `meta.yml`.

## [Unreleased]

Redesign into per-step CPU/GPU modules — see [docs/architecture.md](docs/architecture.md). The README has a table of
the renamed parameters.

### Added

- Samplesheet input (`--input`: `id,rna,atac,metadata`) checked with nf-schema 3.1.0, and parameter validation;
  `--help`, `--help <param>` and `--helpFull`.
- Samples can have RNA, ATAC or both; `--membership_assay` picks the assay the methods run on in paired samples.
- Per-sample cell metadata (samplesheet `metadata` column); its `sample_id` column is optional.
- `memberships.csv` and `ingested.csv` output indexes and `pipeline_info/versions.csv`.
- Profiles `local`, `docker`, `test`, `test_atac`, `test_paired` and `test_10x`.
- Modules `cellgeni/adata/convert` and `cellgeni/adata/attachobs` (adata-cli 0.17.6), with nf-tests.
- `--allow_dtype_change` to accept digit-only categorical obs columns becoming integers when attaching
  metadata (cellgeni/adata-cli#124); without it the run stops.
- `CLAUDE.md` maintainer notes and `docs/architecture.md` (design contracts and roadmap).
- Synthetic test data and its generator (`tests/data/make_test_data.py`): paired RNA+ATAC sample, a second RNA
  sample with overlapping barcodes, counts in `layers/counts`, a 10x multiome `.h5`, a 10x mtx folder and cell
  metadata.
- nf-test configuration and a CI workflow that lints the pipeline with Nextflow 26.10.7.
- GPU job settings `--gpu_steps`, `--gpuqueue`, `--gpu_gmem` and `--costcode` for the GPU-capable steps of the
  coming per-step modules (no step has a GPU path yet).
- Dockerfiles and lock files for the images of the per-step modules (`docker/`): scanpy 1.12.6, rapids-singlecell
  0.18.1 (CUDA 12), snapATAC2 2.10.0 and SEACells a1173c8 (CPU and GPU), all on Python 3.12; GPU smoke tests for
  the farm (`tests/gpu/`).
- `scio.py`, the input/output helper of the step modules (raw directories, slim and attach objects), with tests in
  CI.
- Representation modules with CPU and GPU paths, tested against the v1 in-memory preprocessing (scanpy, muon,
  snapATAC2): `scanpy/normalizetotal`, `scanpy/log1p`, `scanpy/hvg`, `scanpy/pca` (with an in-memory `--scale`
  instead of a separate scale step, which would write a dense layer), `scanpy/neighbors`, `atac/tfidf`, `atac/lsi`,
  `snapatac2/selectfeatures` (CPU) and `snapatac2/spectral` (CuPy port of the cosine embedding on GPU), and
  `anndata/assemble` (a base object plus step outputs in one copy).
- `--atac_embedding lsi|spectral`: ATAC methods can run on snapATAC2's spectral embedding.
- `--rna_counts_key`, `--atac_counts_key` (counts in `X`, `layers/<name>` or `raw/X`), representation options
  (`--normalize_target_sum`, `--hvg_flavor`, `--hvg_batch_key`, `--pca_scale`, `--tfidf_scale_factor`,
  `--lsi_drop_first`, `--select_*`, `--spectral_*`, `--seed`) and step output formats (`--output_mode`,
  `--adata_format`, `--sparse_format`).
- Representation step outputs under `<sample>/<assay>/` and their index `representations.csv`.

### Changed

- Flat parameters: `--methods seacells,hmetacells` instead of `--seacells.enabled`/`--hierarchial.enabled`,
  `--metacell_size` (was `--seacells.gamma`), `--seacells_<option>`, `--hmetacells_<option>`, shared
  `--hvg_n_top_genes`/`--pca_n_comps`/`--lsi_n_comps`/`--embedding_key`, `--constrain_by` and `--eval_labels`
  (were `--celltype_label`), `--barcode_delimiter` (was `--delimiter`).
- 10x inputs (`.h5` files, mtx folders) are detected and converted with adata-cli, replacing
  `ToH5ad`/`AttachCellMetadata`; ATAC keeps its peaks. A bare `.mtx` file is no longer accepted; give its folder.
- Barcode suffixes are applied once, at ingestion (raw inputs used to be suffixed twice).
- The membership assay's embedding is computed by the representation modules, one step per task, and SEACells and
  hmetacells read it (`--precomputed`) instead of preprocessing internally; on the test data hmetacells gives the
  same metacells as 0.1.0.
- Outputs go to `<outdir>/<sample>/<method>/` (Nextflow output block, copied by default; `--publish_mode`), run
  reports to `<outdir>/pipeline_info/`.
- Typed Nextflow (records, typed parameters, strict configuration); requires Nextflow 26.10.6 or later.
- `.gitignore` is now tracked, with general patterns (Nextflow runtime, local results, logs, `.dev/` scratch).

### Removed

- `--filelist`, `--type`, `--raw`, `--celltype_label`, `--delimiter` and the nested `--seacells.*` and
  `--hierarchial.*` parameters.
- `bin/convert_to_h5ad.py`, `bin/attach_annotation.py`, `modules/utils` and `configs/` (now `conf/`).

## [0.1.0] - 2026-11-26

Last release of the legacy pipeline; baseline for regression checks during the redesign.

- `SEACELLS` (`cellgeni/seacells` 1.1.0 on the registry: tie jitter for zero-distance neighbours,
  waypoint-eigenvector clamp, soft assignments for fewer than five metacells, NaN purity for unlabelled
  metacells) and `HMETACELLS` (`cellgeni/hmetacells` 1.0.0), each with in-script preprocessing.
- `ToH5ad` and `AttachCellMetadata` ingestion scripts.
- Module `meta.yml` output names changed from globs to identifiers to pass the current module/v1 schema.
