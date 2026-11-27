"""snapATAC2 step scripts (selectfeatures, spectral on CPU) against snapATAC2 called directly. Needs the snapatac2
environment (docker/snapatac2/requirements.txt); skipped elsewhere."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

snap = pytest.importorskip("snapatac2")
import anndata as ad  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/data/atac.h5ad"
sys.path.insert(0, str(ROOT / "shared/python"))
import scio  # noqa: E402


def run(name, *args):
    path = ROOT / f"modules/cellgeni/snapatac2/{name}/resources/usr/bin/snapatac2_{name}.py"
    spec = importlib.util.spec_from_file_location(f"snapatac2_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main([str(a) for a in args]) == 0


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_selectfeatures_and_spectral_match_snapatac2(workdir):
    run("selectfeatures", "--adata", DATA, "--prefix", "S1", "--n-features", "2000")
    run("spectral", "--adata", DATA, "--prefix", "S1", "--selected", "S1.selectfeatures.raw", "--n-comps", "20")
    got = scio.Reader(DATA, ["S1.selectfeatures.raw", "S1.spectral.raw"])

    reference = ad.read_h5ad(DATA)
    reference.X = reference.X.astype(np.float64)
    snap.pp.select_features(reference, n_features=2000, verbose=False)
    snap.tl.spectral(reference, n_comps=20, random_state=0)
    assert (got.read("var/selected").to_numpy() == reference.var["selected"].to_numpy()).all()
    np.testing.assert_allclose(got.read("uns/spectral_eigenvalue"), reference.uns["spectral_eigenvalue"], rtol=1e-6)
    np.testing.assert_allclose(got.read("obsm/X_spectral"), reference.obsm["X_spectral"], rtol=1e-4, atol=1e-5)


def test_spectral_on_every_feature_with_nystrom(workdir):
    run("spectral", "--adata", DATA, "--prefix", "S1", "--all-features", "--n-comps", "10", "--sample-size", "200")
    embedding = scio.Reader(DATA, ["S1.spectral.raw"]).read("obsm/X_spectral")
    assert embedding.shape[0] == 400 and 1 <= embedding.shape[1] <= 10


def test_gpu_path_refuses_what_it_does_not_port(workdir):
    with pytest.raises(SystemExit, match="full cosine embedding only"):
        run("spectral", "--adata", DATA, "--prefix", "S1", "--all-features", "--sample-size", "100", "--gpu")
