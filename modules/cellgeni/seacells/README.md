# cellgeni/seacells

SEACells metacells from a precomputed embedding: a membership table and a JSON summary. The unified SEACells model
(commit a1173c8) runs on CPU and, when the task has an accelerator, on the device.

```nextflow
nextflow.enable.types = true
include { SEACELLS } from 'cellgeni/seacells'

SEACELLS(channel.of(record(id: 'S1.rna', adata: file('S1_rna.h5ad'), embedding: file('S1.pca.raw'),
                           embedding_key: 'X_pca', output_prefix: null)))
```

## Input record

| Field | Type | Description |
|---|---|---|
| `id` | String | Item id; default output prefix |
| `adata` | Path | Base AnnData object (`.h5ad` or `.zarr`), staged under `input/`; supplies the barcodes |
| `embedding` | Path? | Step output holding `obsm/<embedding_key>` (default: the base object) |
| `embedding_key` | String? | obsm key of the embedding (default `X_pca`) |
| `output_prefix` | String? | Output prefix (default `id`) |

The step output must hold the same cells as the base object, in the same order.

## Output record

`record(id, membership, summary)`:

- `<prefix>.membership.csv` — `barcode,metacell`, one row per cell; an empty `metacell` means unassigned.
- `<prefix>.seacells.json` — `converged`, `n_cells`, `n_metacells`, `requested_metacells`, `metacell_size`
  (min/median/max), `rss`, `embedding_key`, `backend`, `gpu`, `gpu_fallback`, `jitter`, `seed`, `max_iter`
  and tool versions. `gpu` is the device that produced the partition, not the one asked for: a cuSPARSE
  fallback sets `gpu: false` and `gpu_fallback: true`.

## Options

- `ext.args`: one of `--n-metacells <int>` or `--gamma <int>` (cells per metacell) is **required**. Also
  `--n-waypoint-eigs` (default 10), `--n-neighbors` (default 15),
  `--convergence-epsilon` (default 1e-5), `--min-iter` (10), `--max-iter` (50), `--use-sparse` and
  `--seed` (default 0).
- `--use-sparse` selects the legacy sparse CPU backend; SEACells refuses to combine it with the unified model, so
  it cannot be given together with an accelerator.
- GPU: give the task an accelerator (`accelerator = 1`; on LSF also request the GPU, e.g. with clusterOptions).
  The GPU path needs an A100 or H100 (CUDA 12). A cuSPARSE failure falls back to the CPU once, and the
  summary then records the CPU as the device, so a parity check cannot read that result as a GPU one.

## Behaviour worth knowing

- **Exact ties are jittered.** Cells with identical embedding coordinates give zero neighbour distances, which
  sparse matrices do not store; SEACells then builds a zero kernel bandwidth and collapses into one metacell, and
  Palantir indexes past the end of a row. A deterministic jitter is applied only when a preflight graph shows too
  few positive distances, and its scale is recorded in the summary.
- **Non-convergence is reported, not raised.** SEACells sets the labels and only then complains, so the partition
  is usable; `converged: false` in the summary is the signal. With a seed a retry reproduces the same partition,
  so retrying non-convergence cannot help — raise `--max-iter` instead.
- **Randomness is global.** SEACells draws from `np.random`, so `--seed` fixes a run on one backend. CPU and GPU
  give different but equally good partitions at scale: compare RSS, purity and metacell sizes, not ARI.

Containers: `quay.io/cellgeni/seacells:a1173c8-cpu`, `quay.io/cellgeni/seacells:a1173c8-gpu`. Versions are emitted
on the `versions` topic.
