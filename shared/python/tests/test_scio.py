"""Tests for scio: every output mode x object format x sparse format x element type, alignment checks, errors."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import scio  # noqa: E402

N_OBS, N_VARS = 40, 30


def make_adata(n_obs=N_OBS, n_vars=N_VARS, seed=0, reverse_obs=False):
    rng = np.random.default_rng(seed)
    counts = sp.random(n_obs, n_vars, density=0.3, format="csr", dtype=np.float32, rng=rng)
    counts.data = np.ceil(counts.data * 10)
    obs_names = [f"c{i}" for i in range(n_obs)]
    obs = pd.DataFrame({"celltype": pd.Categorical(rng.choice(["a", "b"], n_obs))}, index=obs_names)
    var = pd.DataFrame({"gene_ids": [f"ENSG{i}" for i in range(n_vars)]}, index=[f"g{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=counts, obs=obs, var=var, obsm={"X_old": rng.normal(size=(n_obs, 2))}, uns={"note": "base"})
    adata.layers["counts"] = counts.copy()
    return adata[::-1].copy() if reverse_obs else adata


def write_adata(adata, path):
    adata.write_zarr(path) if path.suffix == ".zarr" else adata.write_h5ad(path)
    return path


def tree_digest(path):
    files = sorted(p for p in Path(path).rglob("*") if p.is_file()) if Path(path).is_dir() else [Path(path)]
    h = hashlib.sha256()
    for f in files:
        h.update(f.read_bytes())
    return h.hexdigest()


def new_elements(n_obs=N_OBS, n_vars=N_VARS):
    rng = np.random.default_rng(1)
    rank = rng.random(n_vars)
    rank[::4] = np.nan
    cluster = pd.Categorical(rng.choice(["x", "y", "z"], n_obs), categories=["z", "y", "x"])
    cluster[0] = np.nan
    return {
        "layers/log1p": sp.random(n_obs, n_vars, density=0.3, format="csr", dtype=np.float32, rng=rng),
        "layers/scaled": rng.normal(size=(n_obs, n_vars)).astype(np.float32),
        "obsm/X_pca": rng.normal(size=(n_obs, 5)).astype(np.float32),
        "varm/PCs": rng.normal(size=(n_vars, 5)),
        "obsp/connectivities": sp.random(n_obs, n_obs, density=0.1, format="csr", dtype=np.float32, rng=rng),
        "obsp/distances": sp.random(n_obs, n_obs, density=0.1, format="csc", dtype=np.float64, rng=rng),
        "obs/n_counts": rng.integers(0, 100, n_obs),
        "obs/cluster": cluster,
        "obs/flag": rng.random(n_obs) > 0.5,
        "obs/label": np.array([f"label {i}" for i in range(n_obs)], dtype=object),
        "var/highly_variable": rng.random(n_vars) > 0.5,
        "var/means": rng.random(n_vars).astype(np.float32),
        "var/rank": rank,
        "var/leiden": pd.Categorical(rng.integers(0, 3, n_vars), categories=[0, 1, 2]),
        "uns/pca": {
            "variance": rng.random(5).astype(np.float32),
            "params": {"n_comps": 5, "zero_center": True, "mask_var": None, "method": "arpack"},
            "labels": np.array(["a", "b"], dtype=object),
        },
    }


def assert_same(key, got, want):
    group, _ = scio.split_key(key)
    if group in scio.FRAME_GROUPS:
        assert isinstance(got, pd.Series)
        pd.testing.assert_series_equal(got.reset_index(drop=True), pd.Series(want), check_names=False)
    elif sp.issparse(want):
        assert sp.issparse(got) and got.format == want.format and got.dtype == want.dtype
        np.testing.assert_array_equal(got.toarray(), want.toarray())
    elif isinstance(want, dict):
        assert set(got) == set(want), key
        for k in want:
            assert_same(key, got[k], want[k])
    elif isinstance(want, np.ndarray) and want.dtype == object:
        assert list(got) == list(want)
    elif isinstance(want, np.ndarray):
        if group != "uns":  # uns arrays go through JSON in raw mode (float32 -> float64)
            assert got.dtype == want.dtype
        np.testing.assert_allclose(got, want, rtol=1e-6)
    else:
        assert got == want


@pytest.fixture(params=["h5ad", "zarr"])
def base(tmp_path, request):
    return write_adata(make_adata(), tmp_path / f"base.{request.param}")


@pytest.fixture
def out_dir(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.chdir(out)
    return out


def write_step(base, step="test", sources=(), **kwargs):
    reader = scio.Reader(base, sources)
    writer = scio.Writer(reader, step, "S1.rna", **kwargs)
    for key, value in new_elements().items():
        writer.add(key, value)
    return writer.write(params={"n_comps": 5, "flavor": None}, versions=scio.versions("anndata"))


def manifest_of(path):
    if path.is_dir() and (path / "manifest.json").exists():
        return json.loads((path / "manifest.json").read_text())
    adata = ad.read_zarr(path) if path.suffix == ".zarr" else ad.read_h5ad(path)
    return {step: json.loads(text) for step, text in adata.uns[scio.UNS_KEY].items()}


# The sparse format applies to raw outputs, the object format to slim and attach outputs
OUTPUTS = [("raw", "h5ad", "npz"), ("raw", "h5ad", "mtx")] + [
    (mode, adata_format, "npz") for mode in ("slim", "attach") for adata_format in scio.ADATA_FORMATS
]


@pytest.mark.parametrize(("mode", "adata_format", "sparse_format"), OUTPUTS)
def test_round_trip(base, out_dir, mode, adata_format, sparse_format):
    out = write_step(base, mode=mode, adata_format=adata_format, sparse_format=sparse_format)
    assert out.name == f"S1.rna.test.{'raw' if mode == 'raw' else adata_format}"
    reader = scio.Reader(base, [out])
    for key, value in new_elements().items():
        assert reader.where(key) == out
        assert_same(key, reader.read(key), value)
    manifest = manifest_of(out)
    manifest = manifest if mode == "raw" else manifest["test"]
    assert manifest["schema"] == scio.SCHEMA and manifest["step"] == "test" and manifest["mode"] == mode
    assert (manifest["n_obs"], manifest["n_vars"]) == (N_OBS, N_VARS)
    assert manifest["obs_names_sha256"] == scio.names_digest(f"c{i}" for i in range(N_OBS))
    assert manifest["params"] == {"n_comps": 5, "flavor": None}
    assert manifest["versions"]["scio"] == scio.__version__ and "anndata" in manifest["versions"]
    assert set(manifest["elements"]) == set(new_elements())
    assert manifest["elements"]["layers/log1p"]["encoding"] == "csr"
    assert manifest["elements"]["var/leiden"]["categories"] == [0, 1, 2]


def test_raw_layout(base, out_dir):
    out = write_step(base, sparse_format="mtx")
    files = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file())
    assert files == [
        "layers/log1p.mtx", "layers/scaled.npy", "manifest.json", "obs.csv", "obsm/X_pca.npy",
        "obsp/connectivities.mtx", "obsp/distances.mtx", "uns/pca.json", "var.csv", "varm/PCs.npy",
    ]
    var = pd.read_csv(out / "var.csv", index_col=0)
    assert var.index.name == "_index" and list(var.index[:2]) == ["g0", "g1"]


def test_attach_keeps_the_base(base, out_dir):
    before = tree_digest(base)
    out = write_step(base, mode="attach", adata_format=base.suffix[1:])
    assert tree_digest(base) == before
    adata = ad.read_zarr(out) if out.suffix == ".zarr" else ad.read_h5ad(out)
    original = make_adata()
    np.testing.assert_array_equal(adata.X.toarray(), original.X.toarray())
    assert {"counts", "log1p", "scaled"} <= set(adata.layers)
    assert {"X_old", "X_pca"} <= set(adata.obsm)
    assert {"celltype", "n_counts", "cluster"} <= set(adata.obs.columns)
    assert list(adata.var["gene_ids"][:2]) == ["ENSG0", "ENSG1"]
    assert adata.uns["note"] == "base" and "pca" in adata.uns


def test_attach_chain_keeps_every_manifest(base, out_dir):
    first = write_step(base, step="first", mode="attach", adata_format="h5ad")
    second = write_step(first, step="second", mode="attach", adata_format="h5ad")
    assert set(manifest_of(second)) == {"first", "second"}


def test_slim_holds_only_the_new_elements(base, out_dir):
    out = write_step(base, mode="slim", adata_format="zarr")
    adata = ad.read_zarr(out)
    assert adata.X is None and adata.shape == (N_OBS, N_VARS)
    assert set(adata.layers) == {"log1p", "scaled"} and set(adata.obsm) == {"X_pca"}
    assert set(adata.obs.columns) == {"n_counts", "cluster", "flag", "label"}


def test_sources_win_over_the_base(base, out_dir):
    reader = scio.Reader(base)
    raw = scio.Writer(reader, "redo", "S1.rna").add("obsm/X_old", np.zeros((N_OBS, 2))).write()
    assert reader.read("obsm/X_old").any()
    assert not scio.Reader(base, [raw]).read("obsm/X_old").any()
    assert scio.Reader(base, [raw]).where("layers/counts") == base


def test_missing_sources_are_skipped(base):
    assert scio.Reader(base, [None, ""]).sources == []


@pytest.mark.parametrize("mode", ["raw", "slim"])
def test_misaligned_sources_are_rejected(base, out_dir, tmp_path, mode):
    out = write_step(base, mode=mode)
    reversed_base = write_adata(make_adata(reverse_obs=True), tmp_path / "reversed.h5ad")
    with pytest.raises(scio.AlignmentError, match="cell names differ"):
        scio.Reader(reversed_base, [out])
    smaller = write_adata(make_adata()[:30].copy(), tmp_path / "smaller.h5ad")
    with pytest.raises(scio.AlignmentError, match="holds 40 cells x 30 features"):
        scio.Reader(smaller, [out])


def test_raw_counts_must_match_var(tmp_path):
    adata = make_adata()
    adata.raw = adata
    good = write_adata(adata, tmp_path / "good.h5ad")
    assert scio.Reader(good).read("raw/X").shape == (N_OBS, N_VARS)
    adata = make_adata()
    adata.raw = adata
    bad = write_adata(adata[:, :20].copy(), tmp_path / "bad.h5ad")
    with pytest.raises(ValueError, match="same features as var"):
        scio.Reader(bad).read("raw/X")


def test_store_keys(base, out_dir):
    store = scio.open_store(base)
    assert {"X", "layers/counts", "obsm/X_old", "obs/celltype", "var/gene_ids", "uns/note"} <= set(store.keys())
    assert set(scio.open_store(write_step(base)).keys()) == set(new_elements())


def test_missing_key_lists_where_it_looked(base):
    with pytest.raises(KeyError, match="obsm/X_nothing is in none of"):
        scio.Reader(base).read("obsm/X_nothing")


def test_writer_checks_what_it_is_given(base, out_dir):
    writer = scio.Writer(scio.Reader(base), "check", "S1.rna")
    with pytest.raises(ValueError, match="expected 40 x any"):
        writer.add("obsm/X_pca", np.zeros((39, 2)))
    with pytest.raises(ValueError, match="two-dimensional"):
        writer.add("obsm/X_pca", np.zeros(40))
    with pytest.raises(ValueError, match="expected 40 x 40"):
        writer.add("obsp/connectivities", sp.csr_matrix((40, 30)))
    with pytest.raises(ValueError, match="never write X"):
        writer.add("X", np.zeros((40, 30)))
    with pytest.raises(ValueError, match="unsupported element key"):
        writer.add("obsm/a/b", np.zeros((40, 2)))
    with pytest.raises(TypeError, match="cannot store a DataFrame"):
        writer.add("uns/frame", pd.DataFrame({"a": [1]}))
    with pytest.raises(TypeError, match="not supported"):
        writer.add("obs/when", pd.to_datetime(np.arange(40), unit="D"))
    writer.add("obsm/X_pca", np.zeros((40, 2)))
    with pytest.raises(ValueError, match="added twice"):
        writer.add("obsm/X_pca", np.zeros((40, 2)))
    with pytest.raises(ValueError, match="step name"):
        scio.Writer(scio.Reader(base), "Bad-Step", "S1")
    with pytest.raises(ValueError, match="not one of"):
        scio.Writer(scio.Reader(base), "check", "S1", mode="copy")
    with pytest.raises(ValueError, match="nothing to write"):
        scio.Writer(scio.Reader(base), "empty", "S1").write()


def test_columns_are_positional(base, out_dir):
    values = pd.Series(np.arange(N_OBS), index=[f"other{i}" for i in reversed(range(N_OBS))])
    out = scio.Writer(scio.Reader(base), "pos", "S1").add("obs/rank", values).write()
    got = scio.Reader(base, [out]).read("obs/rank")
    assert list(got.index[:2]) == ["c0", "c1"] and list(got[:3]) == [0, 1, 2]


def test_open_store_errors(tmp_path):
    (tmp_path / "S1.pca.raw").mkdir()
    with pytest.raises(ValueError, match="without manifest.json"):
        scio.open_store(tmp_path / "S1.pca.raw")
    (tmp_path / "table.csv").write_text("a\n")
    with pytest.raises(ValueError, match="expected a raw step directory"):
        scio.open_store(tmp_path / "table.csv")
    (tmp_path / "other.raw").mkdir()
    (tmp_path / "other.raw" / "manifest.json").write_text('{"schema": "something else"}')
    with pytest.raises(ValueError, match="not an nf-metacells step manifest"):
        scio.open_store(tmp_path / "other.raw")


def test_check_nnz():
    matrix = sp.random(10, 10, density=0.5, format="csr", rng=0)
    scio.check_nnz(matrix)
    with pytest.raises(ValueError, match="run this step on CPU"):
        scio.check_nnz(matrix, limit=10)


def test_names_digest_is_sha256_of_newline_joined_names():
    assert scio.names_digest(["a", "b"]) == hashlib.sha256(b"a\nb").hexdigest()


def test_io_arguments():
    parser = scio.add_io_arguments(argparse.ArgumentParser())
    args = parser.parse_args(["--adata", "x.h5ad", "--prefix", "S1.rna"])
    assert (args.output_mode, args.adata_format, args.sparse_format, args.gpu) == ("raw", "h5ad", "npz", False)
