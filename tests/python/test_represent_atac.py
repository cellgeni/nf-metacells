"""ATAC TF-IDF and LSI step scripts on CPU against muon 0.1.7's tfidf and lsi, which the v1 modules ran."""

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy.sparse import csr_matrix, dia_matrix
from scipy.sparse.linalg import svds

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/data/atac.h5ad"
sys.path.insert(0, str(ROOT / "shared/python"))
import scio  # noqa: E402


def run(tool, name, *args):
    path = ROOT / f"modules/cellgeni/{tool}/{name}/resources/usr/bin/{tool}_{name}.py"
    spec = importlib.util.spec_from_file_location(f"{tool}_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main([str(a) for a in args]) == 0


def muon_tfidf(counts, scale_factor=1e4):
    """muon 0.1.7 atac.pp.tfidf with its defaults (log_tf, log_idf), as the v1 modules called it."""
    n_peaks = np.asarray(counts.sum(axis=1)).reshape(-1)
    tf = np.dot(dia_matrix((1.0 / n_peaks, 0), shape=(n_peaks.size, n_peaks.size)), counts)
    tf = np.log1p(tf * scale_factor)
    idf = np.log1p(np.asarray(counts.shape[0] / counts.sum(axis=0)).reshape(-1))
    return np.nan_to_num(np.dot(tf, dia_matrix((idf, 0), shape=(idf.size, idf.size))), nan=0.0)


def muon_lsi(x, n_comps):
    """muon 0.1.7 atac.tl.lsi (scale_embeddings=True)."""
    u, s, vt = svds(x, k=n_comps)
    u, s, vt = u[:, ::-1], s[::-1], vt[::-1, :]
    return (u - u.mean(axis=0)) / u.std(axis=0), s / np.sqrt(x.shape[0] - 1), vt.T


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_tfidf_and_lsi_match_muon(workdir):
    run("atac", "tfidf", "--adata", DATA, "--prefix", "S1")
    run("atac", "lsi", "--adata", DATA, "--prefix", "S1", "--tfidf", "S1.tfidf.raw", "--n-comps", "30")
    got = scio.Reader(DATA, ["S1.tfidf.raw", "S1.lsi.raw"])

    # v1 removed peaks without counts before muon; ours keeps them as zero columns
    counts = ad.read_h5ad(DATA).X.astype(np.float64)
    keep = np.asarray(counts.sum(axis=0)).ravel() > 0
    reference = csr_matrix(muon_tfidf(counts[:, keep]))
    tfidf = got.read("layers/tfidf")
    assert tfidf.dtype == np.float32 and tfidf[:, ~keep].nnz == 0
    np.testing.assert_allclose(tfidf[:, keep].toarray(), reference.toarray(), rtol=1e-5, atol=1e-6)

    embedding, stdev, loadings = muon_lsi(reference, 30)
    lsi = got.read("obsm/X_lsi")
    for i in range(30):  # the sign of each component is arbitrary in muon
        assert abs(np.corrcoef(lsi[:, i], embedding[:, i])[0, 1]) > 0.9999
    np.testing.assert_allclose(got.read("uns/lsi")["stdev"], stdev, rtol=1e-4)
    assert got.read("varm/LSI").shape == (counts.shape[1], 30)
    assert not got.read("varm/LSI")[~keep].any()


def test_lsi_options(workdir):
    run("atac", "tfidf", "--adata", DATA, "--prefix", "S1", "--output-mode", "slim")
    run("atac", "lsi", "--adata", DATA, "--prefix", "S1", "--tfidf", "S1.tfidf.h5ad", "--n-comps", "10",
        "--drop-first", "--no-scale-embeddings")
    got = scio.Reader(DATA, ["S1.lsi.raw"])
    lsi = got.read("obsm/X_lsi")
    assert lsi.shape == (400, 10)
    np.testing.assert_allclose(np.linalg.norm(lsi, axis=0), 1, rtol=1e-4)  # raw singular vectors
    # signs follow svd_flip: each component's largest entry is positive
    assert (lsi[np.abs(lsi).argmax(axis=0), np.arange(10)] > 0).all()


def test_log_tfidf_needs_plain_terms(workdir):
    with pytest.raises(SystemExit, match="needs --no-log-tf"):
        run("atac", "tfidf", "--adata", DATA, "--prefix", "S1", "--log-tfidf")
    run("atac", "tfidf", "--adata", DATA, "--prefix", "S1", "--log-tfidf", "--no-log-tf", "--no-log-idf")
