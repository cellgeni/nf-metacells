"""CPU/GPU parity of the RNA step scripts (normalizetotal, log1p, hvg, pca, neighbors) in the rapids-singlecell image.

Runs each step script on CPU and with --gpu on synthetic counts (20000 cells x 3000 genes, 8 cell types), chained
through raw step directories as the pipeline chains them, and checks the plan's thresholds: layers within 1e-5,
HVG Jaccard >= 0.95, PCA explained-variance ratios within 1e-3 (relative) with |correlation| >= 0.99 for the
leading components, kNN overlap >= 0.95 (GPU graph vs an exact CPU search). Prints one SMOKE_REPORT JSON line;
exits non-zero when a check fails.

    python tests/gpu/parity_rna.py        (with shared/python on PYTHONPATH)
"""

import importlib.util
import json
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
import scio  # noqa: E402  (shared/python on PYTHONPATH)

report = {"check": "parity_rna", "ok": True}


def step(name):
    path = ROOT / f"modules/cellgeni/scanpy/{name}/resources/usr/bin/scanpy_{name}.py"
    spec = importlib.util.spec_from_file_location(f"scanpy_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def chain(base, prefix, gpu):
    flag = ["--gpu"] if gpu else []
    timings = {}
    for name, args in (
        ("normalizetotal", ["--target-sum", "1e4"]),
        ("log1p", ["--normalized", f"{prefix}.normalizetotal.raw"]),
        ("hvg", ["--log1p", f"{prefix}.log1p.raw", "--n-top-genes", "2000"]),
        ("pca", ["--matrix", f"{prefix}.log1p.raw", "--mask", f"{prefix}.hvg.raw", "--n-comps", "50"]),
        ("neighbors", ["--embedding", f"{prefix}.pca.raw"]),
    ):
        t = time.time()
        assert step(name).main(["--adata", str(base), "--prefix", prefix, *args, *flag]) == 0
        timings[name] = round(time.time() - t, 2)
    outputs = [f"{prefix}.{s}.raw" for s in ("normalizetotal", "log1p", "hvg", "pca", "neighbors")]
    return scio.Reader(base, outputs), timings


def neighbour_sets(reader):
    """Each cell's neighbours in obsp/distances, without the cell itself."""
    d = reader.read("obsp/distances").tocsr()
    return [set(d.indices[d.indptr[i] : d.indptr[i + 1]]) - {i} for i in range(d.shape[0])]


def exact_sets(embedding, k):
    """Each cell's k nearest neighbours by exact Euclidean search, without the cell itself."""
    from sklearn.neighbors import NearestNeighbors

    _, idx = NearestNeighbors(n_neighbors=k + 1, algorithm="brute").fit(embedding).kneighbors(embedding)
    return [set(row[row != i][:k]) for i, row in enumerate(idx)]


def overlap(sets, reference):
    """Mean share of each cell's neighbours that are in its reference neighbours (truncated to the same size)."""
    return round(float(np.mean([len(a & b) / max(len(a), 1) for a, b in zip(sets, reference)])), 4)


def check(name, value, ok):
    report[name] = value
    report["ok"] = report["ok"] and bool(ok)


with tempfile.TemporaryDirectory() as tmp:
    import os

    os.chdir(tmp)
    rng = np.random.default_rng(0)
    n_cells, n_genes, n_types = 20000, 3000, 8
    types = rng.integers(0, n_types, n_cells)
    means = rng.gamma(0.6, 1.5, size=(n_types, n_genes))
    counts = sp.csr_matrix(rng.poisson(means[types] * rng.uniform(0.5, 2, n_cells)[:, None]).astype(np.int32))
    names = {"obs": [f"c{i}" for i in range(n_cells)], "var": [f"g{i}" for i in range(n_genes)]}
    ad.AnnData(counts, obs=dict(index=names["obs"]), var=dict(index=names["var"])).write_h5ad("base.h5ad")

    cpu, report["cpu_s"] = chain("base.h5ad", "cpu", gpu=False)
    gpu, report["gpu_s"] = chain("base.h5ad", "gpu", gpu=True)

    for layer in ("layers/normalized", "layers/log1p"):
        a, b = cpu.read(layer), gpu.read(layer)
        diff = float(abs(a - b).max() / max(abs(a).max(), 1e-12))
        check(f"{layer}_max_rel_diff", diff, diff <= 1e-5)

    hv_cpu, hv_gpu = (r.read("var/highly_variable").to_numpy(bool) for r in (cpu, gpu))
    jaccard = float((hv_cpu & hv_gpu).sum() / (hv_cpu | hv_gpu).sum())
    check("hvg_jaccard", round(jaccard, 4), jaccard >= 0.95)

    ratio_cpu, ratio_gpu = (np.asarray(r.read("uns/pca")["variance_ratio"]) for r in (cpu, gpu))
    rel = float(np.max(np.abs(ratio_cpu - ratio_gpu) / ratio_cpu))
    check("pca_variance_ratio_max_rel_diff", rel, rel <= 1e-3)
    x_cpu, x_gpu = cpu.read("obsm/X_pca"), gpu.read("obsm/X_pca")
    # Components with nearly equal eigenvalues can rotate into each other: compare the leading 10
    corr = [abs(np.corrcoef(x_cpu[:, i], x_gpu[:, i])[0, 1]) for i in range(10)]
    check("pca_min_abs_corr_top10", round(float(min(corr)), 5), min(corr) >= 0.99)

    # The neighbors step alone, on the CPU embedding: the GPU graph (exact search) against an exact CPU search.
    # scanpy's own CPU graph is approximate above a few thousand cells (pynndescent; on these data about a third of
    # the cells have their k-th and (k+1)-th neighbours within 0.1% in distance), so its recall and the end-to-end
    # overlap (each path on its own PCA) are reported for information only.
    assert step("neighbors").main(["--adata", "base.h5ad", "--prefix", "gpuknn", "--embedding", "cpu.pca.raw", "--gpu"]) == 0
    embedding = np.asarray(cpu.read("obsm/X_pca"), dtype=np.float32)
    gpu_sets = neighbour_sets(scio.Reader("base.h5ad", ["gpuknn.neighbors.raw"]))
    exact = exact_sets(embedding, max(len(g) for g in gpu_sets))
    check("knn_overlap", overlap(gpu_sets, exact), overlap(gpu_sets, exact) >= 0.95)
    cpu_sets = neighbour_sets(cpu)
    report["knn_recall_cpu_approximate"] = overlap(cpu_sets, exact)
    report["knn_overlap_end_to_end"] = overlap(cpu_sets, neighbour_sets(gpu))

print("SMOKE_REPORT " + json.dumps(report), flush=True)
sys.exit(0 if report["ok"] else 1)
