"""RNA representation step scripts on CPU, chained through raw step directories as the pipeline chains them.

The chain normalizetotal -> log1p -> hvg -> pca -> neighbors must reproduce scanpy's in-memory preprocessing, which is
what the v1 SEACells and hmetacells modules ran internally (normalize_total 1e4, log1p, seurat HVGs, PCA on them).
"""

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/data"
sys.path.insert(0, str(ROOT / "shared/python"))
import scio  # noqa: E402


def step(name):
    """The step script of modules/cellgeni/scanpy/<name> as a module (its scio copy equals shared/python/scio.py)."""
    path = ROOT / f"modules/cellgeni/scanpy/{name}/resources/usr/bin/scanpy_{name}.py"
    spec = importlib.util.spec_from_file_location(f"scanpy_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(name, *args):
    assert step(name).main([str(a) for a in args]) == 0


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def rna_chain(adata, *extra_normalize):
    run("normalizetotal", "--adata", adata, "--prefix", "S1", "--target-sum", "1e4", *extra_normalize)
    run("log1p", "--adata", adata, "--prefix", "S1", "--normalized", "S1.normalizetotal.raw")
    run("hvg", "--adata", adata, "--prefix", "S1", "--log1p", "S1.log1p.raw", "--n-top-genes", "500")
    run("pca", "--adata", adata, "--prefix", "S1", "--matrix", "S1.log1p.raw", "--mask", "S1.hvg.raw", "--n-comps", "30")
    run("neighbors", "--adata", adata, "--prefix", "S1", "--embedding", "S1.pca.raw")
    outputs = [f"S1.{s}.raw" for s in ("normalizetotal", "log1p", "hvg", "pca", "neighbors")]
    return scio.Reader(adata, outputs)


def reference():
    adata = ad.read_h5ad(DATA / "rna.h5ad")
    adata.X = adata.X.astype(np.float32)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=500, flavor="seurat")
    sc.pp.pca(adata, n_comps=30, mask_var="highly_variable", random_state=0)
    sc.pp.neighbors(adata, random_state=0)
    return adata


def test_chain_matches_in_memory_scanpy(workdir):
    got = rna_chain(DATA / "rna.h5ad")
    ref = reference()
    np.testing.assert_allclose(got.read("layers/log1p").toarray(), ref.X.toarray(), rtol=1e-6)
    assert (got.read("var/highly_variable").to_numpy() == ref.var["highly_variable"].to_numpy()).all()
    np.testing.assert_allclose(got.read("var/dispersions_norm").to_numpy(), ref.var["dispersions_norm"].to_numpy(), rtol=1e-5)
    np.testing.assert_allclose(got.read("obsm/X_pca"), ref.obsm["X_pca"], rtol=1e-4, atol=1e-4)
    np.testing.assert_allclose(got.read("varm/PCs"), ref.varm["PCs"], rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(got.read("uns/pca")["variance_ratio"], ref.uns["pca"]["variance_ratio"], rtol=1e-5)
    np.testing.assert_allclose(
        got.read("obsp/connectivities").toarray(), ref.obsp["connectivities"].toarray(), rtol=1e-4, atol=1e-6
    )
    assert got.read("uns/neighbors")["params"]["n_neighbors"] == 15


def test_counts_from_a_layer_give_the_same_result(workdir):
    from_x = rna_chain(DATA / "rna.h5ad").read("obsm/X_pca")
    (workdir / "layer").mkdir()
    import os

    os.chdir(workdir / "layer")
    counts = ad.read_h5ad(DATA / "rna_layercounts.h5ad")
    assert "counts" in counts.layers
    from_layer = rna_chain(DATA / "rna_layercounts.h5ad", "--counts-key", "layers/counts").read("obsm/X_pca")
    np.testing.assert_allclose(from_layer, from_x, rtol=1e-5, atol=1e-5)


def test_attach_mode_chain_keeps_everything(workdir):
    base = DATA / "rna.h5ad"
    run("normalizetotal", "--adata", base, "--prefix", "S1", "--output-mode", "attach", "--adata-format", "zarr")
    run("log1p", "--adata", "S1.normalizetotal.zarr", "--prefix", "S1", "--output-mode", "attach", "--adata-format", "zarr")
    adata = ad.read_zarr("S1.log1p.zarr")
    assert {"normalized", "log1p"} <= set(adata.layers)
    assert sp.issparse(adata.X) and adata.X.sum() == ad.read_h5ad(base).X.sum()
    assert set(adata.uns["nf_metacells"]) == {"normalizetotal", "log1p"}


def test_seurat_v3_reads_the_counts(workdir):
    run("hvg", "--adata", DATA / "rna.h5ad", "--prefix", "S1", "--flavor", "seurat_v3", "--n-top-genes", "300")
    got = scio.Reader(DATA / "rna.h5ad", ["S1.hvg.raw"])
    assert got.read("var/highly_variable").sum() == 300
    assert "var/variances_norm" in scio.open_store("S1.hvg.raw").keys()


def test_pca_options(workdir):
    base = DATA / "rna.h5ad"
    run("normalizetotal", "--adata", base, "--prefix", "S1")
    run("log1p", "--adata", base, "--prefix", "S1", "--normalized", "S1.normalizetotal.raw")
    run("pca", "--adata", base, "--prefix", "S1", "--matrix", "S1.log1p.raw", "--no-mask", "--scale",
        "--scale-max-value", "10", "--n-comps", "10000")
    got = scio.Reader(base, ["S1.pca.raw"])
    manifest = scio.open_store("S1.pca.raw").manifest
    n_comps = min(got.n_obs, got.n_vars) - 1
    assert got.read("obsm/X_pca").shape == (got.n_obs, n_comps)
    assert manifest["params"]["n_comps"] == n_comps and manifest["params"]["scale"] is True


def test_too_few_genes_stop_hvg(workdir):
    with pytest.raises(SystemExit, match="highly variable genes"):
        step("hvg").main(["--adata", str(DATA / "rna.h5ad"), "--prefix", "S1", "--flavor", "seurat_v3",
                          "--n-top-genes", "1"])
