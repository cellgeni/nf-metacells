# Changelog

All notable changes to nf-metacells. Pipeline releases use calendar versions (`YY-DDD`) from the redesign
onwards; modules under `modules/cellgeni/` use semver and have their own versions in `meta.yml`.

## [Unreleased]

Redesign into per-step CPU/GPU modules — see [docs/architecture.md](docs/architecture.md).

### Added

- `CLAUDE.md` maintainer notes and `docs/architecture.md` (design contracts and roadmap).
- Synthetic test data and its generator (`tests/data/make_test_data.py`): paired RNA+ATAC sample, a second RNA
  sample with overlapping barcodes, counts in `layers/counts`, a 10x multiome `.h5`, a 10x mtx folder and cell
  metadata.
- nf-test configuration and a CI workflow that lints the pipeline with Nextflow 26.10.7.

- Modules `cellgeni/adata/convert` and `cellgeni/adata/attachobs` (adata-cli 0.17.6), with nf-tests.
- `--allow_dtype_change` to accept digit-only categorical obs columns becoming integers when attaching
  metadata (cellgeni/adata-cli#124); without it the run stops.

### Changed

- Ingestion uses adata-cli instead of `ToH5ad`/`AttachCellMetadata`. `--raw --type atac` keeps peaks.
- `--delimiter` is applied once at ingestion for every input (raw inputs used to be suffixed twice).
- Filelist paths may be relative to the launch directory.
- `.gitignore` is now tracked, with general patterns (Nextflow runtime, local results, logs, `.dev/` scratch).

### Removed

- `bin/convert_to_h5ad.py`, `bin/attach_annotation.py`, `modules/utils`. A bare `.mtx` file is no longer accepted
  as raw input; give its 10x folder.

## [0.1.0] - 2026-11-26

Last release of the legacy pipeline; baseline for regression checks during the redesign.

- `SEACELLS` (`cellgeni/seacells` 1.1.0 on the registry: tie jitter for zero-distance neighbours,
  waypoint-eigenvector clamp, soft assignments for fewer than five metacells, NaN purity for unlabelled
  metacells) and `HMETACELLS` (`cellgeni/hmetacells` 1.0.0), each with in-script preprocessing.
- `ToH5ad` and `AttachCellMetadata` ingestion scripts.
- Module `meta.yml` output names changed from globs to identifiers to pass the current module/v1 schema.
