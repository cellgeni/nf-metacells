# nf-metacells — maintainer notes

Nextflow pipeline that builds metacells (SEACells, hierarchical; more methods planned) from scRNA-seq and
scATAC-seq AnnData. It is being redesigned into per-step CPU/GPU modules: the contracts and the roadmap are in
[docs/architecture.md](docs/architecture.md). Tag `0.1.0` is the last release of the legacy pipeline and the
regression baseline.

## Environment (Sanger farm)

```bash
module load cellgen/nextflow/26.10.7 cellgen/nf-test cellgen/uv
```

- GPU jobs: LSF queue `gpu` (A100-80GB, H100-80GB, V100); driver 550 = CUDA 12.4. rapids-singlecell ≥ 0.15
  wheels have no sm_70 kernels, so GPU requests ask for `gmem` ≥ 40000 MB, which only A100/H100 satisfy.
- GitHub CLI: `/software/cellgen/cellgeni/github-cli/2.93.0/bin/gh` (not on PATH).
- Registry: `nextflow module validate|publish modules/cellgeni/<name>` (needs `NXF_REGISTRY_TOKEN`);
  `nextflow module spec` needs `-namespace cellgeni`.
- Scratch work (spikes, GPU smoke jobs, baselines) goes in `.dev/` — git-ignored and on lustre, so GPU nodes see it.

## Commands

- Lint: `nextflow lint .` — must report 0 errors.
- Module tests: `nf-test test modules/cellgeni/<tool>/<step> --profile singularity`.
- Test data: `tests/data/make_test_data.py` (its header has the `uv run` command). Synthetic only — never
  commit real data.

## Conventions

- Every new Nextflow script uses strict syntax and `nextflow.enable.types = true`. A typed process takes one
  record and emits one record; records are declared inline in modules (included types cannot be aliased),
  pipeline record types live in `workflows/types.nf`.
- Every channel item has a unique `id` (`<sample>.<assay>[.<method>.<run>]`); joins across items use `sample`.
- One output field per processing step. In raw mode a step writes `<id>.<step>.raw/`, laid out like the AnnData
  store (`obsm/X_pca.npy`, `uns/pca.json`, …) plus `manifest.json`; slim/attach modes write one `.h5ad|.zarr`.
- GPU: modules branch on `task.accelerator` (container + `--gpu` flag) and never set it; the pipeline config
  sets `accelerator` from `--gpu_steps` and translates it to LSF `-gpu` options and Singularity `--nv`.
- Containers: pinned `quay.io/cellgeni/*` tags, never `latest`; Dockerfiles live in `docker/<image>/`.
- Versions: `topic:` records, `record(process: task.process, tool: '…', version: eval('…')) >> 'versions'`.
- Params are flat (`seacells_gamma`, not `seacells.gamma`) and kept in sync between `params {}` and
  `nextflow_schema.json`.
- meta.yml input/output names must be identifiers; the module/v1 schema rejects glob names like `'*.csv'`.

## Known gotchas

- `join(by: 'id')` warns "Join fields … overwritten" when a step fills a field the state already has (as
  null). That is the intended right-wins update; don't "fix" it with `subMap`, which loses static types.
- Untyped workflows can call typed processes only by passing `record(...)` values — maps fail at runtime.
- `nextflow module spec` regenerates glob output names that `nextflow module validate` then rejects.
- adata-cli 0.17.6: `import dataframe` turns digit-only categoricals (e.g. leiden) into int64 even with
  `--categorical` (cellgeni/adata-cli#124); bare `.mtx` files are refused (10x folders work); no `.npz`
  import/export; never use `--inplace` on Nextflow-staged inputs (it follows the symlink).
- SEACells: pin commit a1173c8 (unified GPU model, `use_unified=True`); `use_sparse` is CPU-only; its `gpu` extra
  pins CUDA 13 packages, so images install the cu12 RAPIDS stack instead.
- cuSPARSE/cupyx use int32 indices: matrices with ≥ 2^31 non-zeros silently truncate on GPU — check before
  moving data to the device.
- pandas 3 makes anndata ≥ 0.12 write string indices/columns as `nullable-string-array`, which older anndata
  (e.g. the legacy 0.9.2 images) and other readers can't open. Images and the test-data generator pin
  `pandas<3` (snapATAC2 2.10 requires it anyway) until the ecosystem catches up.
- `cellgen/nextflow/26.10.7` is the small launcher, which downloads the framework into `~/.nextflow` on first use.

## Git and releases

- One branch per phase or feature, merged into `dev`; `main` only moves at releases (fast-forward).
- Pipeline releases use calendar versions (`YY-DDD`, as nf-cellbender); modules use semver and are published
  from `modules/cellgeni/` to the `cellgeni` registry namespace. Published versions are never changed — bump.
- Record every user-visible change in `CHANGELOG.md`.
