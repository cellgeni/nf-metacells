# quay.io/cellgeni/rapids-singlecell:0.18.1: GPU scverse stack (CUDA 12, RAPIDS 26.6) for nf-metacells
# (see docker/README.md). The CUDA libraries come from the wheels; the driver comes from the host.
ARG PYTHON_IMAGE=python:3.12.16-slim-bookworm

FROM ${PYTHON_IMAGE} AS build
RUN pip install --no-cache-dir uv==0.13.5
COPY requirements.txt /tmp/requirements.txt
RUN uv venv /env --python /usr/local/bin/python3.12 \
    && uv pip install --python /env/bin/python --no-cache --compile-bytecode -r /tmp/requirements.txt

FROM ${PYTHON_IMAGE}
LABEL org.opencontainers.image.source="https://github.com/cellgeni/nf-metacells" \
      org.opencontainers.image.description="rapids-singlecell 0.18.1, RAPIDS 26.6 cu12, scanpy 1.12.6 (GPU) for nf-metacells"
# procps: Nextflow reads task metrics with ps
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /env /env
# Images are read-only under Singularity: numba, CuPy and matplotlib caches go to /tmp
ENV PATH="/env/bin:$PATH" \
    PYTHONNOUSERSITE=1 \
    NUMBA_CACHE_DIR=/tmp/numba \
    CUPY_CACHE_DIR=/tmp/cupy \
    MPLCONFIGDIR=/tmp/matplotlib \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility
COPY Dockerfile requirements.in requirements.txt /docker/
