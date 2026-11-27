"""GPU smoke test of the seacells GPU image: the unified model fits on the GPU and matches the CPU fit.

Synthetic data: 3000 cells in 8 Gaussian clusters on a 20-dimensional embedding, 40 metacells. The same seed must
give the same partition on CPU and GPU (as in the Phase 0 spike). Prints one SMOKE_REPORT JSON line; exits non-zero
when a check fails.

    python tests/gpu/smoke_seacells.py
"""

import json
import sys
import time
import warnings

import anndata as ad
import numpy as np
from sklearn.metrics import adjusted_rand_score

warnings.filterwarnings("ignore")
import SEACells  # noqa: E402

report = {"image": "seacells-gpu", "ok": True}
rng = np.random.default_rng(0)
n_cells, n_dims, n_clusters, k = 3000, 20, 8, 40
centres = rng.normal(scale=4, size=(n_clusters, n_dims))
clusters = rng.integers(0, n_clusters, n_cells)
embedding = (centres[clusters] + rng.normal(size=(n_cells, n_dims))).astype(np.float32)


def fit(gpu):
    adata = ad.AnnData(np.zeros((n_cells, 1), dtype=np.float32), obsm={"X_pca": embedding})
    np.random.seed(0)  # SEACells draws from the global numpy generator
    t = time.time()
    model = SEACells.core.SEACells(
        adata, build_kernel_on="X_pca", n_SEACells=k, n_waypoint_eigs=10, convergence_epsilon=1e-5,
        use_unified=True, use_gpu=gpu, verbose=False,
    )
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    try:
        model.fit(min_iter=10, max_iter=50)
        converged = True
    except RuntimeWarning:  # raised after the labels are set
        converged = False
    labels = adata.obs["SEACell"].astype(str).to_numpy()
    return labels, {"converged": converged, "seconds": round(time.time() - t, 1), "n_metacells": len(set(labels))}


try:
    import cupy as cp

    props = cp.cuda.runtime.getDeviceProperties(0)
    report.update(device=props["name"].decode(), driver=cp.cuda.runtime.driverGetVersion())
    gpu_labels, report["gpu"] = fit(gpu=True)
    cpu_labels, report["cpu"] = fit(gpu=False)
    report["ari_cpu_gpu"] = round(float(adjusted_rand_score(cpu_labels, gpu_labels)), 4)
    report["purity"] = round(
        float(np.mean([np.bincount(clusters[gpu_labels == m]).max() / (gpu_labels == m).sum() for m in set(gpu_labels)])), 4
    )
    report["ok"] = report["ari_cpu_gpu"] == 1.0 and report["gpu"]["n_metacells"] == k and report["purity"] > 0.95
except Exception as error:
    report.update(ok=False, error=f"{type(error).__name__}: {str(error)[:300]}")
print("SMOKE_REPORT " + json.dumps(report), flush=True)
sys.exit(0 if report["ok"] else 1)
