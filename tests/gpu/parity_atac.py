"""CPU/GPU parity of the ATAC step scripts (atac/tfidf, atac/lsi, snapatac2/spectral) on synthetic peak counts.

Two parts, because the CPU spectral embedding needs snapATAC2 and the GPU path needs the rapids-singlecell image:

    python tests/gpu/parity_atac.py prepare DIR    (snapatac2 image) data + CPU spectral embedding
    python tests/gpu/parity_atac.py compare DIR    (rapids-singlecell image, GPU) TF-IDF and LSI both ways, GPU
                                                   spectral, then the checks; prints one SMOKE_REPORT JSON line

Thresholds (plan): TF-IDF within 1e-5; LSI stdev within 1e-3 (relative) and |correlation| >= 0.99 for the leading
10 components; spectral eigenvalues within 1e-3 (relative) and subspace similarity >= 0.99 over the shared
components. Run with shared/python on PYTHONPATH.
"""

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
import scio  # noqa: E402  (shared/python on PYTHONPATH)


def step(tool, name, *args):
    path = ROOT / f"modules/cellgeni/{tool}/{name}/resources/usr/bin/{tool}_{name}.py"
    spec = importlib.util.spec_from_file_location(f"{tool}_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    started = time.time()
    assert module.main([str(a) for a in args]) == 0
    return round(time.time() - started, 2)


def prepare(workdir):
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(0)
    n_cells, n_peaks, n_types = 20000, 50000, 8
    types = rng.integers(0, n_types, n_cells)
    rates = rng.gamma(0.3, 0.1, size=(n_types, n_peaks))
    depth = rng.uniform(0.5, 2, n_cells)[:, None]
    # In chunks of 1000 cells: the dense rate matrix of all cells would take 8 GB
    counts = sp.vstack([
        sp.csr_matrix(rng.poisson(rates[types[i : i + 1000]] * depth[i : i + 1000]).astype(np.int32))
        for i in range(0, n_cells, 1000)
    ]).tocsr()
    names = [f"chr1:{i * 500 + 1}-{(i + 1) * 500}" for i in range(n_peaks)]
    ad.AnnData(counts, obs=dict(index=[f"c{i}" for i in range(n_cells)]), var=dict(index=names)).write_h5ad("base.h5ad")
    timings = {"selectfeatures": step("snapatac2", "selectfeatures", "--adata", "base.h5ad", "--prefix", "cpu",
                                      "--n-features", "20000", "--threads", "4")}
    timings["spectral"] = step("snapatac2", "spectral", "--adata", "base.h5ad", "--prefix", "cpu", "--selected",
                               "cpu.selectfeatures.raw", "--n-comps", "30", "--threads", "4")
    Path("cpu_timings.json").write_text(json.dumps(timings))
    print(f"prepared {workdir}: {n_cells} cells x {n_peaks} peaks, nnz {counts.nnz}", flush=True)


def compare(workdir):
    report = {"check": "parity_atac", "ok": True, "cpu_s": json.loads(Path("cpu_timings.json").read_text()), "gpu_s": {}}

    def check(name, value, ok):
        report[name] = value
        report["ok"] = report["ok"] and bool(ok)

    for prefix, gpu in (("cpu", []), ("gpu", ["--gpu"])):
        seconds = step("atac", "tfidf", "--adata", "base.h5ad", "--prefix", prefix, *gpu)
        seconds_lsi = step("atac", "lsi", "--adata", "base.h5ad", "--prefix", prefix, "--tfidf", f"{prefix}.tfidf.raw",
                           "--n-comps", "30", *gpu)
        report["cpu_s" if not gpu else "gpu_s"].update(tfidf=seconds, lsi=seconds_lsi)
    report["gpu_s"]["spectral"] = step("snapatac2", "spectral", "--adata", "base.h5ad", "--prefix", "gpu", "--selected",
                                       "cpu.selectfeatures.raw", "--n-comps", "30", "--gpu")
    cpu = scio.Reader("base.h5ad", ["cpu.tfidf.raw", "cpu.lsi.raw", "cpu.spectral.raw"])
    gpu = scio.Reader("base.h5ad", ["gpu.tfidf.raw", "gpu.lsi.raw", "gpu.spectral.raw"])

    a, b = cpu.read("layers/tfidf"), gpu.read("layers/tfidf")
    check("tfidf_max_rel_diff", float(abs(a - b).max() / abs(a).max()), abs(a - b).max() / abs(a).max() <= 1e-5)

    sd_cpu, sd_gpu = np.asarray(cpu.read("uns/lsi")["stdev"]), np.asarray(gpu.read("uns/lsi")["stdev"])
    check("lsi_stdev_max_rel_diff", float(np.max(np.abs(sd_cpu - sd_gpu) / sd_cpu)),
          np.max(np.abs(sd_cpu - sd_gpu) / sd_cpu) <= 1e-3)
    x_cpu, x_gpu = cpu.read("obsm/X_lsi"), gpu.read("obsm/X_lsi")
    corr = [abs(np.corrcoef(x_cpu[:, i], x_gpu[:, i])[0, 1]) for i in range(10)]
    check("lsi_min_abs_corr_top10", round(float(min(corr)), 5), min(corr) >= 0.99)

    ev_cpu, ev_gpu = np.asarray(cpu.read("uns/spectral_eigenvalue")), np.asarray(gpu.read("uns/spectral_eigenvalue"))
    k = min(ev_cpu.size, ev_gpu.size)
    report["spectral_components"] = [int(ev_cpu.size), int(ev_gpu.size)]
    check("spectral_eigenvalue_max_rel_diff", float(np.max(np.abs(ev_cpu[:k] - ev_gpu[:k]) / np.abs(ev_cpu[:k]))),
          np.max(np.abs(ev_cpu[:k] - ev_gpu[:k]) / np.abs(ev_cpu[:k])) <= 1e-3)
    q_cpu, _ = np.linalg.qr(cpu.read("obsm/X_spectral")[:, :k].astype(np.float64))
    q_gpu, _ = np.linalg.qr(gpu.read("obsm/X_spectral")[:, :k].astype(np.float64))
    cosines = np.linalg.svd(q_cpu.T @ q_gpu, compute_uv=False)  # cosines of the principal angles
    check("spectral_subspace_similarity", round(float(cosines.mean()), 5), cosines.mean() >= 0.99)
    print("SMOKE_REPORT " + json.dumps(report), flush=True)
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    command, workdir = sys.argv[1], Path(sys.argv[2])
    workdir.mkdir(parents=True, exist_ok=True)
    os.chdir(workdir)
    sys.exit(prepare(workdir) if command == "prepare" else compare(workdir))
