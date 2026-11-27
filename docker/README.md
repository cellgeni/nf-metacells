# Container images

Modules pin these images on quay.io/cellgeni by immutable tags; none uses `latest`.

| Image | Tag | Contents | Used by |
|---|---|---|---|
| `scanpy` | `1.12.6` | Python 3.12, scanpy 1.12.6, anndata 0.13.4, scikit-misc, scikit-network | CPU path of scanpy/\*, atac/\*, hmetacells, aggregate, evaluate, assemble |
| `rapids-singlecell` | `0.18.1` | the `scanpy` stack plus rapids-singlecell 0.18.1, RAPIDS 26.6 (cuDF, cuML, cuGraph, cuVS) and CuPy 14.2, CUDA 12 | GPU path of the same modules and of snapatac2/spectral |
| `snapatac2` | `2.10.0` | Python 3.12, snapATAC2 2.10.0, anndata 0.13.4 | snapatac2/selectfeatures, CPU path of snapatac2/spectral |
| `seacells` | `a1173c8-cpu`, `a1173c8-gpu` | SEACells at commit a1173c8, palantir 1.4.6, scanpy 1.12.6; `-gpu` adds cuML 26.6 and CuPy 14.2 | seacells |

`Dockerfile_python` and `Dockerfile_seacells` build the `latest` images of the v1 SEACells and hmetacells modules;
they go when those modules are replaced (roadmap phase 5).

## Rules

- Every image runs Python 3.12 with anndata 0.13.4, scanpy 1.12.6 and pandas < 3, so CPU and GPU results compare.
- `requirements*.in` lists what an image needs; `requirements*.txt` is the lock the build installs. The SEACells
  locks are the environments of the Phase 0 spike. After editing an `.in` file, recompile its lock:

  ```bash
  uv pip compile docker/<image>/requirements.in -o docker/<image>/requirements.txt \
      --python-version 3.12 --python-platform x86_64-manylinux_2_28 --no-header
  ```

- GPU images use CUDA 12 builds. The farm's driver (550) supports CUDA 12.4 and runs the CUDA 12.9 wheels through
  minor-version compatibility; CUDA 13 needs driver 580 or newer. The CUDA libraries come from the wheels and the
  driver from the host (Singularity `--nv`, Docker `--gpus all`), so every image starts from `python:3.12-slim`.
- Tags name the tool version (for SEACells, its commit and variant). Once a module uses a tag it is never overwritten:
  a rebuild with the same tool version is pushed as `<tag>-2`, `<tag>-3`, …
- Images are read-only under Singularity, so the numba, CuPy and matplotlib caches go to `/tmp`.

## Build and push (Sanger farm)

The cellgeni Docker daemon builds and pushes images but does not run containers:

```bash
module load cellgen/docker
docker build -t quay.io/cellgeni/scanpy:1.12.6 docker/scanpy
docker build --build-arg VARIANT=gpu -t quay.io/cellgeni/seacells:a1173c8-gpu docker/seacells
docker push quay.io/cellgeni/<image>:<tag>
```

Repositories created on quay.io are private: make them public, so that the farm and CI can pull without logging in.

## Test the GPU images

`tests/gpu/run_smoke.sh` submits one LSF job per GPU model (A100 and H100 by default) that runs
`tests/gpu/smoke_rapids.py` and `tests/gpu/smoke_seacells.py` in the GPU images. Run it after changing a GPU image.
