# Container images

Modules pin every image by an immutable tag; none uses `latest`. The CPU images are Seqera Containers (Wave)
builds of the conda environment files here, with their Python packages installed from pip; the CUDA and SEACells
images are our own builds on quay.io/cellgeni.

| Image | Contents | Used by |
|---|---|---|
| `community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90` ([scanpy](scanpy/environment.yml)) | Python 3.12, scanpy 1.12.6, anndata 0.13.4, pandas 2.3.3, scikit-misc, scikit-network | CPU path of scanpy/\*, atac/\*, hmetacells, aggregate, evaluate, assemble |
| `community.wave.seqera.io/library/python_pip_gcc_conda-gcc-specs_pruned:b8199b0e3c90fd81` ([snapatac2](snapatac2/environment.yml)) | Python 3.12, snapATAC2 2.10.0, anndata 0.13.4, pandas 2.3.3 | snapatac2/selectfeatures, CPU path of snapatac2/spectral |
| `quay.io/cellgeni/rapids-singlecell:0.18.1` | the scanpy stack plus rapids-singlecell 0.18.1, RAPIDS 26.6 (cuDF, cuML, cuGraph, cuVS) and CuPy 14.2, CUDA 12 | GPU path of the same modules and of snapatac2/spectral |
| `quay.io/cellgeni/seacells:a1173c8-cpu`, `a1173c8-gpu` | SEACells at commit a1173c8, palantir 1.4.6, scanpy 1.12.6; `-gpu` adds cuML 26.6 and CuPy 14.2 | seacells |

`Dockerfile_python` and `Dockerfile_seacells` build the `latest` images of the v1 SEACells and hmetacells modules;
they go when those modules are replaced (roadmap phase 5).

## Rules

- Every image runs Python 3.12 with anndata 0.13.4, scanpy 1.12.6 and pandas < 3, so CPU and GPU results compare.
- Wave images set no environment beyond `PATH`, so module scripts export `PYTHONNOUSERSITE=1` (packages in
  `~/.local` must not shadow the image's) and put the numba and matplotlib caches under `$TMPDIR` (images are
  read-only under Singularity).
- GPU images use CUDA 12 builds. The farm's driver (550) supports CUDA 12.4 and runs the CUDA 12.9 wheels through
  minor-version compatibility; CUDA 13 needs driver 580 or newer. The CUDA libraries come from the wheels and the
  driver from the host (Singularity `--nv`, Docker `--gpus all`), so every GPU image starts from `python:3.12-slim`.
- On the farm, Nextflow pulls images into `/nfs/cellgeni/singularity/images`. To pull one by hand, keep
  `SINGULARITY_TMPDIR` on local disk: a temporary directory on lustre gave an unpadded squashfs that the kernel
  refused to mount ("bad superblock").

## Wave images (CPU)

`<image>/environment.yml` pins conda-forge Python 3.12 and pip, and the Python packages under `pip:`. Wave freezes
each build into the public `community.wave.seqera.io` registry, tagged with a hash of the build: any change to the
file, comments included, gives a new tag. Build after editing a file and copy the printed names into the modules,
their READMEs, this table, `.github/workflows/ci.yml` and `tests/gpu/run_smoke.sh`:

```bash
module load cellgen/nextflow/26.10.7
cd docker/wave && nextflow inspect main.nf -concretize   # 5-10 min per new build; unchanged files return at once
```

PyPI publishes macs3, a snapATAC2 dependency, as source only, so the snapatac2 environment adds conda-forge `gcc`,
`conda-gcc-specs` and `zlib` for pip to compile it.

## quay.io images (CUDA, SEACells)

- `requirements*.in` lists what an image needs; `requirements*.txt` is the lock the build installs. The SEACells
  locks are the environments of the Phase 0 spike. After editing an `.in` file, recompile its lock:

  ```bash
  uv pip compile docker/<image>/requirements.in -o docker/<image>/requirements.txt \
      --python-version 3.12 --python-platform x86_64-manylinux_2_28 --no-header
  ```

- Tags name the tool version (for SEACells, its commit and variant). Once a module uses a tag it is never overwritten:
  a rebuild with the same tool version is pushed as `<tag>-2`, `<tag>-3`, …
- Images are read-only under Singularity, so their Dockerfiles send the numba, CuPy and matplotlib caches to `/tmp`.

The cellgeni Docker daemon builds and pushes images but does not run containers:

```bash
module load cellgen/docker
docker build --build-arg VARIANT=gpu -t quay.io/cellgeni/seacells:a1173c8-gpu docker/seacells
docker push quay.io/cellgeni/<image>:<tag>
```

Repositories created on quay.io are private: make them public, so that the farm and CI can pull without logging in.

## Test the images

- CPU: CI runs `pytest shared/python tests/python` inside the scanpy image and `tests/python/test_snapatac2.py`
  inside the snapatac2 image, and the module nf-tests with Docker.
- GPU: `tests/gpu/run_smoke.sh` submits one LSF job per GPU model (A100 and H100 by default) that runs
  `tests/gpu/smoke_rapids.py` and `tests/gpu/smoke_seacells.py` in the GPU images, and the CPU/GPU parity checks.
  Run it after changing a GPU image.
