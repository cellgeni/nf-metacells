"""anndata/assemble: base object + step outputs (raw directories and objects) -> one object with every element."""

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/data/rna.h5ad"
sys.path.insert(0, str(ROOT / "shared/python"))
import scio  # noqa: E402


def run(tool, name, *args):
    path = ROOT / f"modules/cellgeni/{tool}/{name}/resources/usr/bin/{tool}_{name}.py"
    spec = importlib.util.spec_from_file_location(f"{tool}_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main([str(a) for a in args]) == 0


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.mark.parametrize("adata_format", ["h5ad", "zarr"])
def test_assemble_raw_and_slim_outputs(workdir, adata_format):
    run("scanpy", "hvg", "--adata", DATA, "--prefix", "S1", "--flavor", "seurat_v3", "--n-top-genes", "300",
        "--output-mode", "slim")
    run("scanpy", "pca", "--adata", DATA, "--prefix", "S1", "--matrix-key", "X", "--mask", "S1.hvg.h5ad",
        "--n-comps", "10")
    run("anndata", "assemble", "--adata", DATA, "--prefix", "S1", "--sources", "S1.hvg.h5ad", "S1.pca.raw",
        "--adata-format", adata_format)
    out = f"S1.assemble.{adata_format}"
    adata = ad.read_zarr(out) if adata_format == "zarr" else ad.read_h5ad(out)
    base = ad.read_h5ad(DATA)
    assert adata.shape == base.shape and (adata.X != base.X).nnz == 0
    assert list(adata.obs.columns) == list(base.obs.columns)
    assert adata.var["highly_variable"].sum() == 300 and "variances_norm" in adata.var
    np.testing.assert_allclose(adata.obsm["X_pca"], scio.Reader(DATA, ["S1.pca.raw"]).read("obsm/X_pca"))
    assert set(adata.uns["nf_metacells"]) == {"hvg", "pca", "assemble"}


def test_assemble_without_sources_copies_the_base(workdir):
    run("anndata", "assemble", "--adata", DATA, "--prefix", "S1")
    adata = ad.read_h5ad("S1.assemble.h5ad")
    assert adata.shape == (400, 1000) and set(adata.uns["nf_metacells"]) == {"assemble"}
