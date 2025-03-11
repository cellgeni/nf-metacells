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
