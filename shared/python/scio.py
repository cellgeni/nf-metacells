#!/usr/bin/env python3
"""Shared input/output for the nf-metacells step modules.

Edit shared/python/scio.py only: scripts/sync_shared.sh copies it into each module's resources/usr/bin/ (modules are
published on their own, so each carries a copy) and CI fails when a copy drifts.

Element keys follow the AnnData store: layers/<name>, obsm/<name>, varm/<name>, obsp/<name>, varp/<name>,
obs/<column>, var/<column>, uns/<key>, and X or raw/X for reading counts.

Reading: Reader(base, sources) looks a key up in the sources first (raw step directories, slim or full AnnData
stores, in the given order) and then in the base object. Every source must hold the same cells and features as the
base, in the same order; anything else fails before a step computes on misaligned data.

Writing: Writer(reader, step, prefix, mode) collects the new elements of a step and writes them in one of three modes

    raw     <prefix>.<step>.raw/          one file per element, laid out like the AnnData store, and manifest.json
    slim    <prefix>.<step>.h5ad|.zarr    obs/var names and the new elements only
    attach  <prefix>.<step>.h5ad|.zarr    a copy of the base object plus the new elements

Raw files: dense arrays .npy; sparse matrices .npz (uncompressed) or .mtx; obs/var columns go to obs.csv/var.csv;
uns entries to uns/<key>.json. The manifest records the shape, digests of the obs and var names, every element's
file, encoding and dtype, the step parameters and tool versions. Slim and attach objects keep the same manifest as
a JSON string in uns['nf_metacells'][<step>].

Command line (modules ship scio.py on their PATH):

    scio.py version <distribution>    installed version, for the versions topic of a module
    scio.py describe <path>           elements of a step output or AnnData store, with shapes and dtypes
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import logging
import platform
import re
import shutil
import sys
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator, Iterable

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

__version__ = "1.0.0"

SCHEMA = "nf-metacells/step-output"
SCHEMA_VERSION = 1
MODES = ("raw", "slim", "attach")
ADATA_FORMATS = ("h5ad", "zarr")
SPARSE_FORMATS = ("npz", "mtx")
MATRIX_GROUPS = ("layers", "obsm", "varm", "obsp", "varp")
FRAME_GROUPS = ("obs", "var")
UNS_KEY = "nf_metacells"
# cuSPARSE and cupyx.scipy.sparse index with int32: bigger matrices are silently truncated on the GPU
NNZ_LIMIT = 2**31 - 1

log = logging.getLogger("scio")


class AlignmentError(ValueError):
    """A source holds other cells or features than the base object, or the same ones in another order."""


# Element keys ---------------------------------------------------------------------------------------------------


def split_key(key: str) -> tuple[str, str]:
    """'obsm/X_pca' -> ('obsm', 'X_pca'); 'X' -> ('X', ''); 'raw/X' -> ('raw', 'X')."""
    if key == "X":
        return "X", ""
    group, _, name = key.partition("/")
    if not name or "/" in name or not (group in (*MATRIX_GROUPS, *FRAME_GROUPS, "uns") or key == "raw/X"):
        raise ValueError(
            f"unsupported element key {key!r}: use X, raw/X, layers/<name>, obsm/<name>, varm/<name>, obsp/<name>, "
            "varp/<name>, obs/<column>, var/<column> or uns/<key>"
        )
    return group, name


def _check_shape(key: str, shape: tuple[int, ...], n_obs: int, n_vars: int) -> None:
    group, _ = split_key(key)
    expected = {
        "X": (n_obs, n_vars), "raw": (n_obs, n_vars), "layers": (n_obs, n_vars),
        "obsm": (n_obs, None), "varm": (n_vars, None), "obsp": (n_obs, n_obs), "varp": (n_vars, n_vars),
        "obs": (n_obs,), "var": (n_vars,),
    }[group]
    if len(shape) != len(expected) or any(e is not None and e != s for e, s in zip(expected, shape)):
        want = " x ".join("any" if e is None else str(e) for e in expected)
        hint = ": counts in raw/X must have the same features as var" if group == "raw" else ""
        raise ValueError(f"{key} has shape {tuple(shape)}, expected {want}{hint}")


def names_digest(names: Iterable[Any]) -> str:
    """SHA-256 of the names joined by newlines: what the manifests record for obs and var names."""
    return hashlib.sha256("\n".join(map(str, names)).encode()).hexdigest()


# Stores ---------------------------------------------------------------------------------------------------------


def open_store(path: str | Path) -> RawStore | AnnDataStore:
    """A raw step directory (with manifest.json), an .h5ad file or a .zarr store."""
    path = Path(path)
    if path.is_dir() and (path / "manifest.json").is_file():
        return RawStore(path)
    if path.is_file() and path.suffix == ".h5ad":
        return AnnDataStore(path, "h5ad")
    if path.is_dir() and (path.suffix == ".zarr" or (path / "zarr.json").is_file() or (path / ".zgroup").is_file()):
        return AnnDataStore(path, "zarr")
    if path.is_dir():
        raise ValueError(f"{path} is a directory without manifest.json: raw step outputs need their manifest")
    raise ValueError(f"{path}: expected a raw step directory, an .h5ad file or a .zarr store")


class AnnDataStore:
    """An .h5ad file or .zarr store, read element by element."""

    def __init__(self, path: Path, kind: str):
        self.path = path
        self.kind = kind

    @contextmanager
    def open(self, mode: str = "r") -> Generator[Any, None, None]:
        if self.kind == "h5ad":
            import h5py

            with h5py.File(self.path, mode) as f:
                yield f
        else:
            import zarr

            # Consolidated metadata can be stale in stores edited after writing: read the real metadata
            yield zarr.open_group(self.path, mode=mode, use_consolidated=False)

    def names(self, axis: str) -> np.ndarray:
        with self.open() as f:
            return _read_index(f[axis])

    def has(self, key: str) -> bool:
        group, name = split_key(key)
        with self.open() as f:
            if group == "X":
                return "X" in f
            if group == "raw":
                return "raw" in f and "X" in f["raw"]
            if group not in f:
                return False
            if group in FRAME_GROUPS:
                return name in _columns(f[group])
            return name in f[group]

    def read(self, key: str) -> Any:
        from anndata.io import read_elem

        group, name = split_key(key)
        with self.open() as f:
            if group == "X":
                return read_elem(f["X"])
            if group == "raw":
                return read_elem(f["raw"]["X"])
            return read_elem(f[group][name])

    def keys(self) -> list[str]:
        with self.open() as f:
            keys = ["X"] if "X" in f else []
            for group in (*MATRIX_GROUPS, "uns"):
                if group in f:
                    keys += [f"{group}/{name}" for name in f[group].keys()]
            for group in FRAME_GROUPS:
                keys += [f"{group}/{name}" for name in _columns(f[group])]
            if "raw" in f and "X" in f["raw"]:
                keys.append("raw/X")
        return keys


def _attr(group: Any, name: str) -> Any:
    value = group.attrs.get(name)
    return value.decode() if isinstance(value, bytes) else value


def _read_index(frame_group: Any) -> np.ndarray:
    from anndata.io import read_elem

    return np.asarray(read_elem(frame_group[_attr(frame_group, "_index")]), dtype=object)


def _columns(frame_group: Any) -> list[str]:
    order = _attr(frame_group, "column-order")
    return [] if order is None else [str(c) for c in np.atleast_1d(order)]


class RawStore:
    """A raw step directory: element files and manifest.json."""

    def __init__(self, path: Path):
        self.path = path
        manifest = json.loads((path / "manifest.json").read_text())
        if manifest.get("schema") != SCHEMA:
            raise ValueError(f"{path}/manifest.json is not an nf-metacells step manifest")
        if manifest.get("version", 0) > SCHEMA_VERSION:
            raise ValueError(f"{path} was written by a newer scio (layout version {manifest['version']})")
        self.manifest = manifest
        self._frames: dict[str, pd.DataFrame] = {}

    def has(self, key: str) -> bool:
        return key in self.manifest["elements"]

    def keys(self) -> list[str]:
        return list(self.manifest["elements"])

    def read(self, key: str) -> Any:
        if key not in self.manifest["elements"]:
            raise KeyError(f"{key} is not in {self.path}")
        spec = self.manifest["elements"][key]
        file = self.path / spec["path"]
        encoding = spec["encoding"]
        if encoding == "dense":
            return np.load(file, allow_pickle=False)
        if encoding in ("csr", "csc"):
            matrix = sp.load_npz(file) if file.suffix == ".npz" else scipy.io.mmread(file, spmatrix=False)
            matrix = sp.csr_matrix(matrix) if encoding == "csr" else sp.csc_matrix(matrix)
            return matrix.astype(spec["dtype"], copy=False)
        if encoding == "column":
            if spec["path"] not in self._frames:
                self._frames[spec["path"]] = pd.read_csv(file, dtype=str, keep_default_na=False, index_col=0)
            return _decode_column(self._frames[spec["path"]][split_key(key)[1]], spec)
        if encoding == "json":
            return _from_json(json.loads(file.read_text()))
        raise ValueError(f"{self.path}: unknown encoding {encoding!r} for {key}")


# Columns and JSON ----------------------------------------------------------------------------------------------


def _column_spec(values: pd.Series) -> dict[str, Any]:
    dtype = values.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        return {
            "dtype": "category",
            "categories": _to_json(dtype.categories.to_numpy()),
            "categories_dtype": str(dtype.categories.dtype),
            "ordered": bool(dtype.ordered),
        }
    if (
        pd.api.types.is_bool_dtype(dtype)
        or pd.api.types.is_string_dtype(dtype)
        or (pd.api.types.is_numeric_dtype(dtype) and not pd.api.types.is_complex_dtype(dtype))
    ):
        return {"dtype": str(dtype)}
    raise TypeError(f"column {values.name}: dtype {dtype} is not supported (use numbers, booleans, strings or categories)")


def _decode_column(text: pd.Series, spec: dict[str, Any]) -> Any:
    """CSV strings back to the recorded dtype; 'nan' marks missing values."""
    if spec["dtype"] == "category":
        categories = pd.Index(spec["categories"]).astype(spec["categories_dtype"])
        lookup = {str(c): c for c in categories}
        return pd.Categorical(text.map(lookup), categories=categories, ordered=spec["ordered"])
    dtype = pd.api.types.pandas_dtype(spec["dtype"])
    if pd.api.types.is_bool_dtype(dtype):
        if isinstance(dtype, np.dtype):
            if not text.isin(["True", "False"]).all():
                raise ValueError(f"column {text.name}: missing or invalid values in a bool column")
            return (text == "True").to_numpy()
        return text.map({"True": True, "False": False}).astype(dtype).array
    if pd.api.types.is_numeric_dtype(dtype):
        return pd.Series(pd.to_numeric(text.mask(text == "nan"))).astype(dtype).array
    return text.astype(dtype).array


def _to_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _to_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_json(v) for v in value]
    if isinstance(value, np.ndarray):
        if value.dtype.kind not in "biufUO":
            raise TypeError(f"cannot store an array of dtype {value.dtype} as JSON")
        return _to_json(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    raise TypeError(f"cannot store a {type(value).__name__} in uns: use dicts, lists, arrays, numbers and strings")


def _from_json(value: Any) -> Any:
    """Lists of numbers or strings come back as numpy arrays, as anndata reads uns."""
    if isinstance(value, dict):
        return {k: _from_json(v) for k, v in value.items()}
    if isinstance(value, list):
        try:
            array = np.asarray(value)
        except ValueError:  # ragged
            array = None
        if array is not None and array.dtype.kind in "biuf":
            return array
        if array is not None and array.dtype.kind == "U":
            return array.astype(object)
        return [_from_json(v) for v in value]
    return value


# Reader ---------------------------------------------------------------------------------------------------------


class Reader:
    """The base object and the step outputs a step reads from, checked to hold the same cells and features."""

    def __init__(self, base: str | Path, sources: Iterable[str | Path | None] = ()):
        store = open_store(base)
        if not isinstance(store, AnnDataStore):
            raise ValueError(f"the base object must be an .h5ad file or a .zarr store, not {base}")
        self.base = store
        self.obs_names = store.names("obs")
        self.var_names = store.names("var")
        self.n_obs, self.n_vars = len(self.obs_names), len(self.var_names)
        self.obs_digest = names_digest(self.obs_names)
        self.var_digest = names_digest(self.var_names)
        self.sources: list[RawStore | AnnDataStore] = []
        for source in sources:
            if source is None or str(source) == "":
                continue
            opened = open_store(source)
            self._check_aligned(opened)
            self.sources.append(opened)

    def _check_aligned(self, store: RawStore | AnnDataStore) -> None:
        if isinstance(store, RawStore):
            m = store.manifest
            n_obs, n_vars, obs_digest, var_digest = m["n_obs"], m["n_vars"], m["obs_names_sha256"], m["var_names_sha256"]
        else:
            obs, var = store.names("obs"), store.names("var")
            n_obs, n_vars, obs_digest, var_digest = len(obs), len(var), names_digest(obs), names_digest(var)
        if (n_obs, n_vars) != (self.n_obs, self.n_vars):
            raise AlignmentError(
                f"{store.path} holds {n_obs} cells x {n_vars} features but the base object {self.base.path.name} "
                f"holds {self.n_obs} x {self.n_vars}: it was computed for another object (for example before cells "
                "were filtered)"
            )
        for axis, ours, theirs in (("cell", self.obs_digest, obs_digest), ("feature", self.var_digest, var_digest)):
            if ours != theirs:
                raise AlignmentError(
                    f"{store.path}: its {axis} names differ from the base object {self.base.path.name} "
                    "(same count, other names or another order)"
                )

    def find(self, key: str) -> RawStore | AnnDataStore | None:
        split_key(key)
        for store in (*self.sources, self.base):
            if store.has(key):
                return store
        return None

    def has(self, key: str) -> bool:
        return self.find(key) is not None

    def where(self, key: str) -> Path | None:
        store = self.find(key)
        return None if store is None else store.path

    def read(self, key: str) -> Any:
        """Arrays and sparse matrices as stored; obs/var columns as a Series indexed by the base names."""
        store = self.find(key)
        if store is None:
            searched = ", ".join(str(s.path) for s in (*self.sources, self.base))
            raise KeyError(f"{key} is in none of: {searched}")
        value = store.read(key)
        log.info("read %s from %s", key, store.path)
        group, name = split_key(key)
        if group in FRAME_GROUPS:
            index = self.obs_names if group == "obs" else self.var_names
            values = value.array if isinstance(value, pd.Series) else value
            return pd.Series(values, index=pd.Index(index), name=name)
        if group != "uns":
            _check_shape(key, value.shape, self.n_obs, self.n_vars)
        return value


# Writer ---------------------------------------------------------------------------------------------------------

_STEP_NAME = re.compile(r"^[a-z][a-z0-9_]*$")


class Writer:
    """The new elements of one step, written in raw, slim or attach mode."""

    def __init__(
        self,
        reader: Reader,
        step: str,
        prefix: str,
        mode: str = "raw",
        adata_format: str = "h5ad",
        sparse_format: str = "npz",
    ):
        if not _STEP_NAME.match(step):
            raise ValueError(f"step name {step!r}: use lower-case letters, digits and underscores")
        for value, allowed in ((mode, MODES), (adata_format, ADATA_FORMATS), (sparse_format, SPARSE_FORMATS)):
            if value not in allowed:
                raise ValueError(f"{value!r} is not one of {', '.join(allowed)}")
        self.reader = reader
        self.step = step
        self.prefix = prefix
        self.mode = mode
        self.adata_format = adata_format
        self.sparse_format = sparse_format
        self.elements: dict[str, Any] = {}

    @property
    def path(self) -> Path:
        return Path(f"{self.prefix}.{self.step}.{'raw' if self.mode == 'raw' else self.adata_format}")

    def add(self, key: str, value: Any) -> Writer:
        group, name = split_key(key)
        if group in ("X", "raw"):
            raise ValueError(f"steps never write {key}: write a layer instead")
        if key in self.elements:
            raise ValueError(f"{key} was added twice")
        value = to_cpu(value)
        r = self.reader
        if group in FRAME_GROUPS:
            # Positional: the values follow the base order, whatever index a Series carries
            value = pd.Series(value.array if isinstance(value, pd.Series) else value, name=name)
            _check_shape(key, value.shape, r.n_obs, r.n_vars)
            _column_spec(value)
        elif group == "uns":
            _to_json(value)  # every mode must accept what a step writes, and raw mode stores JSON
        else:
            if sp.issparse(value):
                value = value if value.format in ("csr", "csc") else value.tocsr()
            else:
                value = np.asarray(value)
                if value.dtype.kind not in "biuf":
                    raise TypeError(f"{key}: arrays of dtype {value.dtype} are not supported")
            if group in ("obsm", "varm") and value.ndim != 2:
                raise ValueError(f"{key} must be two-dimensional, got shape {value.shape}")
            _check_shape(key, value.shape, r.n_obs, r.n_vars)
        self.elements[key] = value
        return self

    def write(self, params: dict[str, Any] | None = None, versions: dict[str, str] | None = None) -> Path:
        if not self.elements:
            raise ValueError("nothing to write: add elements first")
        out = self.path
        if out.is_dir():
            shutil.rmtree(out)
        elif out.exists():
            out.unlink()
        r = self.reader
        manifest = {
            "schema": SCHEMA,
            "version": SCHEMA_VERSION,
            "step": self.step,
            "prefix": self.prefix,
            "mode": self.mode,
            "n_obs": r.n_obs,
            "n_vars": r.n_vars,
            "obs_names_sha256": r.obs_digest,
            "var_names_sha256": r.var_digest,
            "elements": {key: _spec(key, value) for key, value in self.elements.items()},
            "params": _to_json(params or {}),
            "versions": dict(versions or {}),
        }
        {"raw": self._write_raw, "slim": self._write_slim, "attach": self._write_attach}[self.mode](out, manifest)
        log.info("wrote %s: %s", out, ", ".join(self.elements))
        return out

    def _write_raw(self, out: Path, manifest: dict[str, Any]) -> None:
        out.mkdir()
        frames: dict[str, dict[str, pd.Series]] = {}
        for key, value in self.elements.items():
            group, name = split_key(key)
            spec = manifest["elements"][key]
            if group in FRAME_GROUPS:
                spec["path"] = f"{group}.csv"
                frames.setdefault(group, {})[name] = value
                continue
            if group == "uns":
                spec["path"] = f"uns/{name}.json"
            elif sp.issparse(value):
                spec["path"] = f"{group}/{name}.{self.sparse_format}"
            else:
                spec["path"] = f"{group}/{name}.npy"
            file = out / spec["path"]
            file.parent.mkdir(parents=True, exist_ok=True)
            if group == "uns":
                file.write_text(json.dumps(_to_json(value), indent=2) + "\n")
            elif not sp.issparse(value):
                np.save(file, value, allow_pickle=False)
            elif self.sparse_format == "npz":
                sp.save_npz(file, value, compressed=False)
            else:
                scipy.io.mmwrite(file, value)
        for group, columns in frames.items():
            names = self.reader.obs_names if group == "obs" else self.reader.var_names
            frame = pd.DataFrame({name: s.array for name, s in columns.items()}, index=pd.Index(names, name="_index"))
            frame.to_csv(out / f"{group}.csv", na_rep="nan")
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    def _write_slim(self, out: Path, manifest: dict[str, Any]) -> None:
        import anndata as ad

        r = self.reader
        adata = ad.AnnData(obs=pd.DataFrame(index=pd.Index(r.obs_names)), var=pd.DataFrame(index=pd.Index(r.var_names)))
        for key, value in self.elements.items():
            group, name = split_key(key)
            if group in FRAME_GROUPS:
                getattr(adata, group)[name] = value.array
            elif group == "uns":
                adata.uns[name] = value
            else:
                getattr(adata, group)[name] = value
        adata.uns[UNS_KEY] = {self.step: json.dumps(manifest)}
        _write_adata(adata, out)

    def _write_attach(self, out: Path, manifest: dict[str, Any]) -> None:
        from anndata.io import read_elem, write_elem

        base = self.reader.base
        if base.kind == self.adata_format and base.kind == "h5ad":
            shutil.copyfile(base.path, out)  # a new file: never edit a (staged, symlinked) input in place
        elif base.kind == self.adata_format:
            shutil.copytree(base.path, out, copy_function=shutil.copyfile)
        else:
            import anndata as ad

            log.warning("attach mode converts %s to %s and loads the whole object", base.path.name, self.adata_format)
            _write_adata(ad.read_h5ad(base.path) if base.kind == "h5ad" else ad.read_zarr(base.path), out)
        frames: dict[str, dict[str, pd.Series]] = {}
        with AnnDataStore(out, self.adata_format).open("r+") as f:
            for key, value in self.elements.items():
                group, name = split_key(key)
                if group in FRAME_GROUPS:
                    frames.setdefault(group, {})[name] = value
                    continue
                parent = _require_group(f, group)
                if name in parent:
                    del parent[name]
                write_elem(parent, name, value)
            for group, columns in frames.items():
                frame = read_elem(f[group])
                for name, values in columns.items():
                    frame[name] = values.array
                del f[group]
                write_elem(f, group, frame)
            uns = _require_group(f, "uns")
            steps = read_elem(uns[UNS_KEY]) if UNS_KEY in uns else {}
            steps[self.step] = json.dumps(manifest)
            if UNS_KEY in uns:
                del uns[UNS_KEY]
            write_elem(uns, UNS_KEY, steps)
        if self.adata_format == "zarr":
            import zarr

            # Readers such as anndata.read_zarr use the consolidated metadata, which the edits above made stale
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Consolidated metadata is currently not part")
                zarr.consolidate_metadata(str(out))


def _spec(key: str, value: Any) -> dict[str, Any]:
    group, _ = split_key(key)
    if group in FRAME_GROUPS:
        return {"encoding": "column", **_column_spec(value)}
    if group == "uns":
        return {"encoding": "json"}
    spec = {"shape": list(value.shape), "dtype": str(value.dtype)}
    if sp.issparse(value):
        return {"encoding": value.format, **spec, "nnz": int(value.nnz)}
    return {"encoding": "dense", **spec}


def _require_group(f: Any, group: str) -> Any:
    if group not in f:
        from anndata.io import write_elem

        write_elem(f, group, {})
    return f[group]


def _write_adata(adata: Any, out: Path) -> None:
    if out.suffix == ".zarr":
        adata.write_zarr(out)
    else:
        adata.write_h5ad(out)


# GPU ------------------------------------------------------------------------------------------------------------


def gpu_init(managed_memory: bool = True) -> None:
    """Route CuPy allocations through RMM; managed memory lets large matrices spill to host memory."""
    import cupy
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    rmm.reinitialize(managed_memory=managed_memory, pool_allocator=False)
    cupy.cuda.set_allocator(rmm_cupy_allocator)


def check_nnz(matrix: Any, limit: int = NNZ_LIMIT) -> None:
    """Refuse sparse matrices whose non-zeros overflow the int32 indices of cuSPARSE."""
    if sp.issparse(matrix) and matrix.nnz > limit:
        raise ValueError(
            f"the matrix has {matrix.nnz} non-zeros, more than cuSPARSE can index ({limit}): run this step on CPU"
        )


def to_gpu(x: Any, dtype: Any = None) -> Any:
    """Copy to the GPU. Integers become float32 unless dtype says otherwise: cupyx sparse has no integer types."""
    import cupy
    import cupyx.scipy.sparse as cpsp

    if sp.issparse(x):
        check_nnz(x)
        x = x.astype(dtype or (np.float32 if x.dtype.kind in "biu" else x.dtype), copy=False)
        return cpsp.csc_matrix(x) if x.format == "csc" else cpsp.csr_matrix(x.tocsr())
    x = np.asarray(x)
    return cupy.asarray(x, dtype=dtype or (np.float32 if x.dtype.kind in "biu" else x.dtype))


def to_cpu(x: Any) -> Any:
    """CuPy arrays and cupyx sparse matrices back to numpy/scipy; anything else unchanged."""
    return x.get() if type(x).__module__.split(".")[0] in ("cupy", "cupyx") else x


# Command line ---------------------------------------------------------------------------------------------------


def add_io_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """The options every step script shares."""
    io = parser.add_argument_group("input and output")
    io.add_argument("--adata", required=True, type=Path, help="base AnnData object (.h5ad or .zarr)")
    io.add_argument("--prefix", required=True, help="output name prefix, usually the record id")
    io.add_argument("--output-mode", choices=MODES, default="raw", help="raw directory, slim object or full copy")
    io.add_argument("--adata-format", choices=ADATA_FORMATS, default="h5ad", help="object format (slim, attach)")
    io.add_argument("--sparse-format", choices=SPARSE_FORMATS, default="npz", help="sparse matrix files (raw)")
    io.add_argument("--gpu", action="store_true", help="compute on the GPU")
    return parser


def writer_from_args(reader: Reader, step: str, args: argparse.Namespace) -> Writer:
    return Writer(reader, step, args.prefix, args.output_mode, args.adata_format, args.sparse_format)


def versions(*packages: str) -> dict[str, str]:
    """Python, scio and the installed versions of the given distributions, for the manifest."""
    found = {"python": platform.python_version(), "scio": __version__}
    for package in packages:
        try:
            found[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    return found


def setup_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(name)s: %(message)s", datefmt="%H:%M:%S")


def describe(path: str | Path) -> list[str]:
    """One line per element of a step output or AnnData store."""
    store = open_store(path)
    lines = []
    if isinstance(store, RawStore):
        m = store.manifest
        lines.append(f"{store.path}: step {m['step']} ({m['mode']}), {m['n_obs']} cells x {m['n_vars']} features")
        for key, spec in m["elements"].items():
            shape = " x ".join(map(str, spec.get("shape", [])))
            lines.append(f"  {key}  {spec['encoding']} {spec.get('dtype', '')} {shape}".rstrip())
        return lines
    obs, var = store.names("obs"), store.names("var")
    lines.append(f"{store.path}: {store.kind}, {len(obs)} cells x {len(var)} features")
    lines += [f"  {key}" for key in store.keys()]
    return lines


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="scio.py", description="nf-metacells step input/output helper")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("version", help="print the installed version of a distribution").add_argument("distribution")
    commands.add_parser("describe", help="list the elements of a step output or AnnData store").add_argument("path")
    args = parser.parse_args(argv)
    if args.command == "version":
        print(importlib.metadata.version(args.distribution))
    else:
        print("\n".join(describe(args.path)))
    return 0


if __name__ == "__main__":
    sys.exit(_main())
