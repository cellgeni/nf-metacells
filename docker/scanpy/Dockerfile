# quay.io/cellgeni/scanpy:1.12.6: CPU scverse stack for nf-metacells (see docker/README.md)
ARG PYTHON_IMAGE=python:3.12.16-slim-bookworm

FROM ${PYTHON_IMAGE} AS build
RUN pip install --no-cache-dir uv==0.13.5
COPY requirements.txt /tmp/requirements.txt
RUN uv venv /env --python /usr/local/bin/python3.12 \
    && uv pip install --python /env/bin/python --no-cache --compile-bytecode -r /tmp/requirements.txt

FROM ${PYTHON_IMAGE}
LABEL org.opencontainers.image.source="https://github.com/cellgeni/nf-metacells" \
      org.opencontainers.image.description="scanpy 1.12.6, anndata 0.13.4 (CPU) for nf-metacells"
# procps: Nextflow reads task metrics with ps
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /env /env
# Images are read-only under Singularity: numba and matplotlib caches go to /tmp
ENV PATH="/env/bin:$PATH" \
    PYTHONNOUSERSITE=1 \
    NUMBA_CACHE_DIR=/tmp/numba \
    MPLCONFIGDIR=/tmp/matplotlib
COPY Dockerfile requirements.in requirements.txt /docker/
