# quay.io/cellgeni/seacells:a1173c8-cpu and :a1173c8-gpu: SEACells at commit a1173c8 for nf-metacells
# (see docker/README.md). VARIANT=gpu adds the cu12 RAPIDS stack (cuML, CuPy) for the unified GPU model.
ARG PYTHON_IMAGE=python:3.12.16-slim-bookworm

FROM ${PYTHON_IMAGE} AS build
ARG VARIANT=cpu
RUN apt-get update && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir uv==0.13.5
COPY requirements-${VARIANT}.txt /tmp/requirements.txt
RUN uv venv /env --python /usr/local/bin/python3.12 \
    && uv pip install --python /env/bin/python --no-cache --compile-bytecode -r /tmp/requirements.txt

FROM ${PYTHON_IMAGE}
ARG VARIANT=cpu
LABEL org.opencontainers.image.source="https://github.com/cellgeni/nf-metacells" \
      org.opencontainers.image.description="SEACells a1173c8 (${VARIANT}) for nf-metacells"
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
COPY Dockerfile requirements-${VARIANT}.in requirements-${VARIANT}.txt /docker/
