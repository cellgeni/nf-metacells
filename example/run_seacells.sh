#!/bin/bash
# Run the SeaCells example

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
