
# nf-metacells

This project uses Nextflow to run SEACells analysis on single-cell RNA/ATAC sequencing data.

## Project Structure

```
/lustre/scratch127/cellgen/cellgeni/aljes/nf-metacells/
├── main.nf
├── nextflow.config
├── example/
│   ├── run_seacells.sh
│   └── samples.csv
└── README.md
```

## Usage

To run the SeaCells analysis, use the provided script in the `example` directory:

```bash
bash example/run_seacells.sh
```

## Parameters

The following parameters can be configured in the `main.nf` script:

- `--filelist`: Path to the CSV file containing sample information.
- `--seacells.enabled`: Enable SeaCells analysis.
- `--seacells.gamma`: Gamma parameter for SeaCells.
- `--seacells.type`: Type of data (e.g., gex).
- `--seacells.n_top_genes`: Number of top genes to consider.
- `--seacells.n_components`: Number of components for dimensionality reduction.
- `--seacells.celltype_label`: Label for cell types.
- `--seacells.convergence_epsilon`: Convergence epsilon for the algorithm.
- `--seacells.min_iterations`: Minimum number of iterations.
- `--seacells.max_iterations`: Maximum number of iterations.

## Example

An example script is provided to demonstrate how to run the analysis:

```bash
#!/bin/bash
# Run the SeaCells example

nextflow run main.nf \
  --filelist samples.csv \
  --seacells.enabled \
  --seacells.gamma 75 \
  --seacells.type gex \
  --seacells.n_top_genes 2000 \
  --seacells.n_components 50 \
  --seacells.celltype_label celltype \
  --seacells.convergence_epsilon 0.00001 \
  --seacells.min_iterations 10 \
  --seacells.max_iterations 50 \
  -resume
```
