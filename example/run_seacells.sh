#!/bin/bash
# SEACells on the samples in samplesheet.csv (id,rna,atac,metadata), with cell-type purity reported

nextflow run cellgeni/nf-metacells -r dev \
  --input samplesheet.csv \
  --methods seacells \
  --metacell_size 75 \
  --hvg_n_top_genes 2000 \
  --pca_n_comps 50 \
  --eval_labels celltype \
  --seacells_convergence_epsilon 0.00001 \
  --seacells_min_iter 10 \
  --seacells_max_iter 50 \
  -resume
