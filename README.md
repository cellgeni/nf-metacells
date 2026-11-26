
# nf-metacells

Pipeline to aggregate metacells/pseudobulks using Nextflow. This pipeline performs metacell-aggregation using SEACells or Hierarchial clustering on single-cell RNA/ATAC sequencing data.

## Usage

```bash
nextflow run main.nf [OPTIONS]
```

### Required Options

- `--filelist` - File containing the list of files to be aggregated
- `--type` - Type of data to be aggregated (gex or atac)
- `--celltype_label` - Cell type label for the aggregation (required for hierarchial clustering; optional for SEACells)

### Optional Options

- `--help` - Show help message
- `--outdir` - Output directory (default `results`)
- `--raw` - Convert 10x `.h5` files or mtx folders to H5AD with [adata-cli](https://github.com/cellgeni/adata-cli). A bare `.mtx` file is not accepted; give its 10x folder
- `--delimiter` - Rename barcodes to `<barcode><delimiter><id>` at ingestion (originals kept in `obs['barcode']`, sample id in `obs['sample']`)
- `--cell_metadata` - CSV with `barcode` and `sample_id` columns whose other columns are attached to `.obs`; cells without metadata are dropped
- `--barcode_column` - `.obs` column matched to the metadata barcodes (default: obs names)
- `--allow_dtype_change` - Accept categorical `.obs` columns with digit-only labels (e.g. leiden) becoming integers when attaching metadata ([adata-cli#124](https://github.com/cellgeni/adata-cli/issues/124)); without it the run stops

### SEACells Options

- `--seacells.enabled` - Enable SEACells aggregation
- `--seacells.n_cells` - Number of metacells to be aggregated
- `--seacells.gamma` - Gamma value for the adaptive bandwidth kernel
- `--seacells.n_top_genes` - Number of top variable genes to be used for aggregation (GEX only). Default: 2000
- `--seacells.n_components` - Number of components to be used for aggregation (PCA for GEX and LSI for ATAC). Default: 50
- `--seacells.convergence_epsilon` - Convergence epsilon for the optimization. Default: 0.00001
- `--seacells.min_iterations` - Minimum number of iterations for the optimization. Default: 10
- `--seacells.max_iterations` - Maximum number of iterations for the optimization. Default: 100
- `--seacells.use_sparse` - Use sparse matrix for the optimization. Default: false
- `--seacells.precomputed` - Use precomputed distance matrix for the aggregation

### Hierarchial Options

- `--hierarchial.enabled` - Enable Hierarchial clustering aggregation
- `--hierarchial.n_min` - Minimum number of clusters for the aggregation
- `--hierarchial.n_max` - Maximum number of clusters for the aggregation
- `--hierarchial.method` - Method for the aggregation (kmeans or louvain)
- `--hierarchial.n_top_genes` - Number of top variable genes to be used for aggregation (GEX only). Default: 2000
- `--hierarchial.n_components` - Number of components to be used for aggregation (PCA for GEX and LSI for ATAC). Default: 50
- `--hierarchial.n_neighbors` - Number of neighbors for the aggregation. Default: 15
- `--hierarchial.precomputed` - Use precomputed distance matrix for the aggregation

## Examples

### 1. Perform metacell aggregation using SEACells

```bash
nextflow run main.nf \
  --filelist filelist.csv \
  --type gex \
  --seacells.enabled \
  --seacells.gamma 75
```

### 2. Full SEACells example with all parameters

```bash
nextflow run main.nf \
  --filelist samples.csv \
  --type gex \
  --seacells.enabled \
  --seacells.gamma 75 \
  --seacells.n_top_genes 2000 \
  --seacells.n_components 50 \
  --celltype_label celltype \
  --seacells.convergence_epsilon 0.00001 \
  --seacells.min_iterations 10 \
  --seacells.max_iterations 50 \
  -resume
```

## Filelist Format

The filelist should be a CSV file with the following format:

```csv
id,path
pbmc_10k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc_10k/pbmc10k/pbmc10k.h5ad
pbmc_3k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc3k/pbmc3k/pbmc3k.h5ad
```

Paths may be absolute or relative to the launch directory.
