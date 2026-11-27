"""CPU/GPU parity of the RNA step scripts (normalizetotal, log1p, hvg, pca, neighbors) in the rapids-singlecell image.

Runs each step script on CPU and with --gpu on synthetic counts (20000 cells x 3000 genes, 8 cell types), chained
through raw step directories as the pipeline chains them, and checks the plan's thresholds: layers within 1e-5,
HVG Jaccard >= 0.95, PCA explained-variance ratios within 1e-3 (relative) with |correlation| >= 0.99 for the
leading components, kNN overlap >= 0.95. Prints one SMOKE_REPORT JSON line; exits non-zero when a check fails.

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

    def neighbour_sets(reader):
        d = reader.read("obsp/distances").tocsr()
        return [set(d.indices[d.indptr[i] : d.indptr[i + 1]]) for i in range(d.shape[0])]

    # The neighbors step alone: GPU graph of the CPU embedding (exact search) vs the CPU graph (approximate above
    # a few thousand cells); end to end (each on its own PCA) is reported for information
    assert step("neighbors").main(["--adata", "base.h5ad", "--prefix", "gpuknn", "--embedding", "cpu.pca.raw", "--gpu"]) == 0
    same = scio.Reader("base.h5ad", ["gpuknn.neighbors.raw"])
    cpu_sets = neighbour_sets(cpu)
    overlap = np.mean([len(a & b) / max(len(a), 1) for a, b in zip(cpu_sets, neighbour_sets(same))])
    check("knn_overlap", round(float(overlap), 4), overlap >= 0.95)
    end_to_end = np.mean([len(a & b) / max(len(a), 1) for a, b in zip(cpu_sets, neighbour_sets(gpu))])
    report["knn_overlap_end_to_end"] = round(float(end_to_end), 4)

print("SMOKE_REPORT " + json.dumps(report), flush=True)
sys.exit(0 if report["ok"] else 1)
