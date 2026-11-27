"""GPU smoke test of the rapids-singlecell image: the CUDA stack runs on the farm driver (550, CUDA 12.4).

Covers the GPU paths nf-metacells uses: rapids-singlecell preprocessing, cuML/cuVS neighbours, a CuPy RawKernel
compiled through NVRTC (as the CuPy spectral port does), cupyx svds/eigsh (LSI, spectral) and the scio GPU helpers.
Prints one SMOKE_REPORT JSON line; exits non-zero when a check fails.

    python tests/gpu/smoke_rapids.py        (with shared/python on PYTHONPATH for the scio checks)
"""

import json
import sys
import time

import anndata as ad
import numpy as np
import scipy.sparse as sp

report = {"image": "rapids-singlecell", "ok": True}
started = time.time()


def check(name, fn):
    t = time.time()
    try:
        value = fn()
        import cupy

        cupy.cuda.Stream.null.synchronize()
        report[name] = round(time.time() - t, 3) if value is None else value
    except Exception as error:  # report every failing check, not only the first
        report[name] = f"FAILED {type(error).__name__}: {str(error)[:200]}"
        report["ok"] = False


rng = np.random.default_rng(0)
n_cells, n_genes, k = 20000, 3000, 8
labels = rng.integers(0, k, n_cells)
counts = sp.csr_matrix(rng.poisson(rng.gamma(0.6, 1.5, size=(k, n_genes))[labels]).astype(np.float32))
adata = ad.AnnData(counts)

import cupy as cp  # noqa: E402
import cupyx.scipy.sparse as cpsp  # noqa: E402
from cupyx.scipy.sparse.linalg import LinearOperator, eigsh, svds  # noqa: E402

props = cp.cuda.runtime.getDeviceProperties(0)
report.update(
    device=props["name"].decode(),
    compute_capability=f"{props['major']}.{props['minor']}",
    driver=cp.cuda.runtime.driverGetVersion(),
    runtime=cp.cuda.runtime.runtimeGetVersion(),
)

try:
    import scio

    def scio_round_trip():
        if (scio.to_cpu(scio.to_gpu(counts)) != counts).nnz:
            raise ValueError("to_gpu/to_cpu changed the matrix")

    check("scio_gpu_init", lambda: scio.gpu_init())
    check("scio_round_trip", scio_round_trip)
except ImportError:
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    rmm.reinitialize(managed_memory=True, pool_allocator=False)
    cp.cuda.set_allocator(rmm_cupy_allocator)
    report["scio"] = "not on PYTHONPATH"

import rapids_singlecell as rsc  # noqa: E402

check("to_gpu", lambda: rsc.get.anndata_to_GPU(adata))
check("normalize_total", lambda: rsc.pp.normalize_total(adata, target_sum=1e4))
check("log1p", lambda: rsc.pp.log1p(adata))
check("hvg_seurat", lambda: rsc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor="seurat"))
check("pca", lambda: rsc.pp.pca(adata, n_comps=30))
check("neighbors", lambda: rsc.pp.neighbors(adata, n_neighbors=15))
report["n_hvg"] = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else None

kernel = cp.RawKernel(
    r"""
extern "C" __global__ void scale_rows(const int* indptr, float* data, const float* s, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) for (int j = indptr[i]; j < indptr[i + 1]; j++) data[j] *= s[i];
}""",
    "scale_rows",
)


def rawkernel():
    m = cpsp.csr_matrix(counts)
    before = float(m.data.sum())
    kernel(((n_cells + 255) // 256,), (256,), (m.indptr.astype(cp.int32), m.data, cp.full(n_cells, 0.5, cp.float32), cp.int32(n_cells)))
    if abs(float(m.data.sum()) - 0.5 * before) > 1e-3 * before:
        raise ValueError("RawKernel result is wrong")


check("rawkernel_nvrtc", rawkernel)
m = cpsp.csr_matrix(counts)
check("svds_k30", lambda: svds(m, k=30) and None)
normalized = m.multiply(1.0 / cp.sqrt(m.multiply(m).sum(axis=1))).tocsr()
operator = LinearOperator((n_cells, n_cells), matvec=lambda v: normalized @ (normalized.T @ v), dtype=cp.float32)
check("eigsh_k30", lambda: eigsh(operator, k=30, which="LM") and None)


def sknetwork():
    import sknetwork  # hmetacells hierarchies (CPU code in the GPU image)

    return sknetwork.__version__


check("sknetwork", sknetwork)
report["total_s"] = round(time.time() - started, 1)
print("SMOKE_REPORT " + json.dumps(report), flush=True)
sys.exit(0 if report["ok"] else 1)
