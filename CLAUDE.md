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

- Lint: `nextflow lint .` — must report 0 errors (it also lints `.dev/`, so delete broken scratch scripts).
- Module tests: `nf-test test modules/cellgeni/<tool>/<step> --profile singularity`.
- Pipeline: `nextflow run . -profile test,local` (also `test_atac`, `test_paired`, `test_10x`; add `-stub-run` for a
  quick wiring check). Test samplesheets use paths relative to the repository root, so run from there.
- Test data: `tests/data/make_test_data.py` (its header has the `uv run` command). Synthetic only — never
  commit real data.
- Shared I/O helper tests: `uv run --no-project --python 3.12 --with-requirements docker/scanpy/requirements.txt
  --with pytest==9.2.0 python -m pytest shared/python -q`.
- GPU images: `tests/gpu/run_smoke.sh` (one LSF job per GPU model); building and tagging rules in `docker/README.md`.

## Conventions

- Every new Nextflow script uses strict syntax and `nextflow.enable.types = true`. A typed process takes one
  record and emits one record; records are declared inline in modules (included types cannot be aliased),
  pipeline record types live in `workflows/types.nf`.
- Every channel item has a unique `id` (`<sample>.<assay>[.<method>.<run>]`); joins across items use `sample`.
- One output field per processing step. In raw mode a step writes `<id>.<step>.raw/`, laid out like the AnnData
  store (`obsm/X_pca.npy`, `uns/pca.json`, …) plus `manifest.json`; slim/attach modes write one `.h5ad|.zarr`.
- GPU: modules branch on `task.accelerator` (container + `--gpu` flag) and never set it. GPU-capable processes
  carry `label 'process_gpu'`; `conf/base.config` sets `accelerator` when `--gpu_steps` names the process' step
  (the last part of its name, lower-cased: `SCANPY_PCA` → `pca`) and adds LSF `-gpu` options and `--nv`.
- Step scripts read and write elements through `scio.py`: edit `shared/python/scio.py` only, then run
  `scripts/sync_shared.sh` (it copies the file into each module's `resources/usr/bin/`; CI checks the copies).
- Containers: pinned `quay.io/cellgeni/*` tags, never `latest`, and no `docker://` prefix (the Docker engine
  rejects it; Singularity adds it itself); Dockerfiles and lock files live in `docker/<image>/`.
- Versions: `topic:` records, `record(process: task.process, tool: '…', version: eval('…')) >> 'versions'`.
- Params are flat (`seacells_min_iter`, not `seacells.min_iter`), typed in `main.nf` `params {}` and documented in
  `nextflow_schema.json`. Only params that config settings use (`outdir`, `publish_mode`) live in `nextflow.config`.
- Layout: `main.nf` (params, entry workflow, output block) → `workflows/metacells.nf` → `subworkflows/local/*`;
  records in `workflows/types.nf`; `conf/base.config` (defaults), `conf/modules.config` (ext.args from params),
  `conf/test*.config`.
- meta.yml input/output names must be identifiers; the module/v1 schema rejects glob names like `'*.csv'`.

## Known gotchas

- `join(by: 'id')` warns "Join fields … overwritten" when a step fills a field the state already has (as
  null). That is the intended right-wins update; don't "fix" it with `subMap`, which loses static types.
- Untyped workflows can call typed processes only by passing `record(...)` values — maps fail at runtime. Typed
  workflows call untyped processes with tuples (`tuple(meta, path)`), which is how the v1 SEACELLS/HMETACELLS run
  until phase 5.
- Typed workflow bodies cannot `return`: use `exit(status, message)`. Scripts cannot have top-level variables:
  use a function. Typed functions read `def name(arg: Type) -> Type { … }`.
- Dynamic process directives need no closure (`container task.accelerator ? GPU_IMAGE : CPU_IMAGE`); `mix` takes
  one channel; a process called with plain values returns a `Value`, which channel operators such as `mix` refuse.
- Config closures (`ext.args`) see the script's params at task time, but in scripts the linter only knows params
  declared in `params {}`.
- nf-schema 3.1.0: its functions take no `workflow` argument (`paramsSummaryLog()`); `samplesheetToList` returns
  rows in the property order of the schema, with empty cells as `[]`; `help` is a `String?` param so that
  `--help <param>` works.
- `nextflow module spec` regenerates glob output names that `nextflow module validate` then rejects.
- adata-cli 0.17.6: `import dataframe` turns digit-only categoricals (e.g. leiden) into int64 even with
  `--categorical` (cellgeni/adata-cli#124); bare `.mtx` files are refused (10x folders work); no `.npz`
  import/export; never use `--inplace` on Nextflow-staged inputs (it follows the symlink). `convert
  --feature-type` fails on 10x input without feature types and on AnnData input (cellgeni/adata-cli#127), so
  `main.nf` guesses typed inputs from the path; a Cell Ranger 2 `.h5` still fails until #127 is fixed.
- SEACells: pin commit a1173c8 (unified GPU model, `use_unified=True`); `use_sparse` is CPU-only; its `gpu` extra
  pins CUDA 13 packages, so images install the cu12 RAPIDS stack instead. It imports IPython without declaring it.
  Identical cells give a zero kernel bandwidth and a silent collapse into one metacell (keep the v1 jitter for
  exact ties); non-convergence raises `RuntimeWarning` after the labels are set; randomness comes from the global
  `np.random`. On GPU, the B update recomputes K·B in every Frank-Wolfe step (~48 s per iteration at 100k cells),
  and cuSPARSE SpGEMM in the kernel step failed once (`CUSPARSE_STATUS_INTERNAL_ERROR`) on a rerun of the same data.
  The GPU is ~9x faster per iteration at 20k cells and ~50x at 100k. With the same seed, CPU and GPU (and the
  dense/sparse backends) give different but equally good partitions at scale, so parity checks compare the
  objective (RSS), purity and metacell sizes, not ARI.
- cuSPARSE/cupyx use int32 indices: matrices with ≥ 2^31 non-zeros silently truncate on GPU — check before
  moving data to the device.
- pandas 3 makes anndata ≥ 0.12 write string indices/columns as `nullable-string-array`, which older anndata
  (e.g. the legacy 0.9.2 images) and other readers can't open. Images and the test-data generator pin
  `pandas<3` (snapATAC2 2.10 requires it anyway) until the ecosystem catches up.
- `cellgen/nextflow/26.10.7` is the small launcher, which downloads the framework into `~/.nextflow` on first use.
  A new shell finds Nextflow 25.04.4 first: load 26.10.7 in the same shell, or `nextflow lint` dies with a
  NullPointerException in `CmdLint.printErrors` instead of listing the errors.

## Git and releases

- One branch per phase or feature, merged into `dev`; `main` only moves at releases (fast-forward).
- Commits are atomic (one concern each) with short, to-the-point one-line messages.
- Pipeline releases use calendar versions (`YY-DDD`, as nf-cellbender); modules use semver and are published
  from `modules/cellgeni/` to the `cellgeni` registry namespace. Published versions are never changed — bump.
- Record every user-visible change in `CHANGELOG.md`.
