
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
- `--output_dir` - Output directory for the aggregated files
- `--raw` - Convert raw files to H5AD format (if raw files' paths are provided)
- `--delimiter` - Specify delimiter if you want to add sample id to the barcode names
- `--cell_metadata` - Specify metadata.csv file to attach to the .obs section of AnnData object
- `--barcode_column` - Barcode column name in the cell metadata file (specify "obs_names" if you want to use the .obs index column)

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
item,filepath
pbmc_10k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc_10k/pbmc10k/pbmc10k.h5ad
pbmc_3k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc3k/pbmc3k/pbmc3k.h5ad
```
