#!/usr/bin/env python3
"""Generate the small synthetic datasets used by nf-metacells tests.

Deterministic (fixed seeds) and structured (cell types with distinct expression and
accessibility programs) so metacell methods have something to find. Regenerate with:

    module load cellgen/uv
    uv run --no-project --python 3.12 --with anndata==0.13.4 --with h5py --with scipy --with "pandas<3" \
        tests/data/make_test_data.py tests/data

pandas is pinned below 3: with pandas 3, anndata >= 0.12 writes string indices as `nullable-string-array`,
which older anndata versions (and the legacy images) cannot read.

Outputs (all small enough to commit):
    rna.h5ad              S1 RNA counts in X (int32), obs: celltype, leiden (digit labels), batch
    atac.h5ad             S1 ATAC peak counts, same barcodes as rna.h5ad (paired), 20 all-zero peaks
    rna_s2.h5ad           S2 RNA counts; shares 50 barcodes with S1 (tests barcode suffixing)
    rna_layercounts.h5ad  S1 RNA with log-normalised X and counts in layers['counts']
    multiome/filtered_feature_bc_matrix.h5   S1 as a 10x multiome .h5 (Gene Expression + Peaks)
    rna_mtx/              S2 RNA as a 10x mtx folder (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
    metadata.csv          barcode,sample_id,annotation,cluster for S1/S2; drops ~10% of S1 cells,
                          adds 3 barcodes absent from the data
    samplesheet_*.csv     inputs for the test profiles (paths relative to the repository root)
"""

import gzip
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

N_TYPES = 5
TYPE_NAMES = np.array(["T", "B", "NK", "Mono", "DC"])
N_GENES = 1000
N_PEAKS = 3000
N_ZERO_PEAKS = 20


def barcodes(rng: np.random.Generator, n: int) -> np.ndarray:
    """Random, unique 10x-style barcodes (16 nt + '-1')."""
    seen: set[str] = set()
    while len(seen) < n:
        seen.add("".join(rng.choice(list("ACGT"), 16)) + "-1")
    return np.array(sorted(seen))


def rna_counts(rng: np.random.Generator, types: np.ndarray) -> sp.csr_matrix:
    """Negative-binomial-like counts with 60 marker genes per cell type."""
    base = rng.gamma(0.3, 0.4, size=N_GENES)
    profiles = np.tile(base, (N_TYPES, 1))
    for t in range(N_TYPES):
        markers = np.arange(t * 60, (t + 1) * 60)
        profiles[t, markers] *= rng.uniform(4, 12, size=markers.size)
    library = rng.lognormal(mean=0.0, sigma=0.35, size=types.size)[:, None]
    mean = profiles[types] * library
    lam = rng.gamma(shape=5.0, scale=mean / 5.0)  # gamma-Poisson mixture
    return sp.csr_matrix(rng.poisson(lam).astype(np.int32))


def atac_counts(rng: np.random.Generator, types: np.ndarray) -> sp.csr_matrix:
    """Sparse accessibility with 200 type-specific peaks per type and empty peaks at the end."""
    rate = np.full((N_TYPES, N_PEAKS), 0.02)
    for t in range(N_TYPES):
        rate[t, t * 200 : (t + 1) * 200] = 0.35
    rate[:, N_PEAKS - N_ZERO_PEAKS :] = 0.0
    depth = rng.lognormal(0.0, 0.3, size=types.size)[:, None]
    return sp.csr_matrix(rng.poisson(rate[types] * depth).astype(np.int32))


def gene_names() -> np.ndarray:
    names = np.array([f"GENE{i:04d}" for i in range(N_GENES)], dtype=object)
    names[[10, 11]] = "DUPGENE"  # duplicated symbol: exercises --make-unique
    return names


def peak_names() -> np.ndarray:
    return np.array([f"chr{1 + i // 1000}:{(i % 1000) * 1000 + 1}-{(i % 1000) * 1000 + 500}" for i in range(N_PEAKS)])


def write_10x_h5(path: Path, obs_names, gene_ids, genes, peaks, rna: sp.csr_matrix, atac: sp.csr_matrix) -> None:
    """Minimal Cell Ranger ARC-style filtered_feature_bc_matrix.h5 (features x cells, CSC)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    mat = sp.vstack([rna.T, atac.T]).tocsc()
    n_feat = len(genes) + len(peaks)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g["barcodes"] = np.array([b.encode() for b in obs_names])
        gz = {"compression": "gzip", "compression_opts": 4}
        g.create_dataset("data", data=mat.data.astype(np.int32), **gz)
        g.create_dataset("indices", data=mat.indices.astype(np.int64), **gz)
        g.create_dataset("indptr", data=mat.indptr.astype(np.int64), **gz)
        g["shape"] = np.array(mat.shape, dtype=np.int32)
        fe = g.create_group("features")
        fe["id"] = np.array([x.encode() for x in list(gene_ids) + list(peaks)])
        fe["name"] = np.array([x.encode() for x in list(genes) + list(peaks)])
        fe["feature_type"] = np.array([b"Gene Expression"] * len(genes) + [b"Peaks"] * len(peaks))
        fe["genome"] = np.array([b"GRCh38"] * n_feat)
        fe["_all_tag_keys"] = np.array([b"genome"])


def write_10x_mtx(folder: Path, obs_names, gene_ids, genes, counts: sp.csr_matrix) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, sp.coo_matrix(counts.T), field="integer")
    with gzip.open(folder / "barcodes.tsv.gz", "wt") as fh:
        fh.write("\n".join(obs_names) + "\n")
    with gzip.open(folder / "features.tsv.gz", "wt") as fh:
        fh.write("".join(f"{i}\t{n}\tGene Expression\n" for i, n in zip(gene_ids, genes)))


def make_obs(obs_names, types: np.ndarray, rng: np.random.Generator) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "celltype": pd.Categorical(TYPE_NAMES[types]),
            "leiden": pd.Categorical(types.astype(str)),  # digit-only labels (cellgeni/adata-cli#124)
            "batch": pd.Categorical(rng.choice(["b1", "b2"], size=types.size)),
        },
        index=pd.Index(obs_names, name=None),
    )


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    genes = gene_names()
    gene_ids = np.array([f"ENSG{i:011d}" for i in range(N_GENES)])
    peaks = peak_names()
    var_rna = pd.DataFrame({"gene_ids": gene_ids}, index=pd.Index(genes, name=None))
    var_rna.index = ad.utils.make_index_unique(var_rna.index.astype(str))
    var_atac = pd.DataFrame(index=pd.Index(peaks, name=None))

    # S1: paired RNA + ATAC, 400 cells
    n1 = 400
    bc1 = barcodes(rng, n1)
    t1 = np.sort(rng.integers(0, N_TYPES, n1))
    rna1, atac1 = rna_counts(rng, t1), atac_counts(rng, t1)
    obs1 = make_obs(bc1, t1, rng)
    ad.AnnData(rna1, obs=obs1, var=var_rna).write_h5ad(out / "rna.h5ad", compression="gzip")
    ad.AnnData(atac1, obs=obs1[["celltype", "leiden"]].copy(), var=var_atac).write_h5ad(out / "atac.h5ad", compression="gzip")

    # S1 with counts kept in layers['counts'] and log-normalised X
    lognorm = rna1.astype(np.float32).multiply(1e4 / np.asarray(rna1.sum(axis=1)).clip(min=1)).tocsr()
    lognorm.data = np.log1p(lognorm.data)
    ad.AnnData(lognorm, obs=obs1, var=var_rna, layers={"counts": rna1}).write_h5ad(
        out / "rna_layercounts.h5ad", compression="gzip"
    )
    write_10x_h5(out / "multiome" / "filtered_feature_bc_matrix.h5", bc1, gene_ids, genes, peaks, rna1, atac1)

    # S2: RNA only, 250 cells, 50 barcodes shared with S1
    n2 = 250
    bc2 = np.array(sorted(set(bc1[:50]) | set(barcodes(np.random.default_rng(1), n2 - 50))))[:n2]
    t2 = np.sort(rng.integers(0, N_TYPES, bc2.size))
    rna2 = rna_counts(rng, t2)
    ad.AnnData(rna2, obs=make_obs(bc2, t2, rng), var=var_rna).write_h5ad(out / "rna_s2.h5ad", compression="gzip")
    write_10x_mtx(out / "rna_mtx", bc2, gene_ids, genes, rna2)

    # Cell metadata: drop ~10% of S1 cells, add 3 unknown barcodes, keep all S2 cells
    keep1 = np.sort(rng.choice(n1, size=int(n1 * 0.9), replace=False))
    meta = pd.concat(
        [
            pd.DataFrame({"barcode": bc1[keep1], "sample_id": "S1", "annotation": TYPE_NAMES[t1[keep1]], "cluster": t1[keep1].astype(str)}),
            pd.DataFrame({"barcode": ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1", "GGGGGGGGGGGGGGGG-1"], "sample_id": "S1", "annotation": "unknown", "cluster": "99"}),
            pd.DataFrame({"barcode": bc2, "sample_id": "S2", "annotation": TYPE_NAMES[t2], "cluster": t2.astype(str)}),
        ],
        ignore_index=True,
    )
    meta.to_csv(out / "metadata.csv", index=False)

    rel = Path("tests/data")
    sheets = {
        "samplesheet_rna.csv": "id,rna\nS1,{d}/rna.h5ad\nS2,{d}/rna_s2.h5ad\n",
        "samplesheet_atac.csv": "id,atac\nS1,{d}/atac.h5ad\n",
        "samplesheet_paired.csv": "id,rna,atac\nS1,{d}/rna.h5ad,{d}/atac.h5ad\n",
        "samplesheet_10x.csv": "id,rna,atac\nS1,{d}/multiome/filtered_feature_bc_matrix.h5,{d}/multiome/filtered_feature_bc_matrix.h5\nS2,{d}/rna_mtx,\n",
    }
    for name, text in sheets.items():
        (out / name).write_text(text.format(d=rel))
    print(f"wrote test data to {out}")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "tests/data"))
