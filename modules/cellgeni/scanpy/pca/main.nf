nextflow.enable.types = true

/*
 * Principal components of a log-normalised layer, on the highly variable genes, into obsm/X_pca, varm/PCs and
 * uns/pca: scanpy on CPU, rapids-singlecell when the task has an accelerator. ext.args: --n-comps <int> (default
 * 50; capped by the data), --matrix-key (default layers/log1p), --mask-key (default var/highly_variable), --no-
 * mask, --no-zero-center, --svd-solver, --scale and --scale-max-value (scale the masked genes in memory first),
 * --seed (default 0).
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process SCANPY_PCA {
    tag id
    label 'process_gpu'
    container task.accelerator ?
        'quay.io/cellgeni/rapids-singlecell:0.18.1' :
        'community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90'

    input:
    record(id: String, adata: Path, matrix: Path?, mask: Path?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, pca: file("${output_prefix ?: id}.pca.*"))

    topic:
    record(
        process: task.process,
        tool: task.accelerator ? 'rapids-singlecell' : 'scanpy',
        version: eval("scio.py version ${task.accelerator ? 'rapids-singlecell-cu12' : 'scanpy'}")
    ) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    export PYTHONNOUSERSITE=1 NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}/numba MPLCONFIGDIR=\${TMPDIR:-/tmp}/matplotlib
    scanpy_pca.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${matrix ? "--matrix ${matrix}" : ''} \\
        ${mask ? "--mask ${mask}" : ''} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.pca.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
