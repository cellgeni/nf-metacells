nextflow.enable.types = true

/*
 * Flags highly variable genes in var/highly_variable, with the flavour's statistics in other var columns; never
 * subsets: scanpy on CPU, rapids-singlecell when the task has an accelerator. ext.args: --flavor
 * seurat|cell_ranger|seurat_v3|seurat_v3_paper (default seurat), --n-top-genes <int>, --batch-key <obs column>,
 * --min-mean, --max-mean, --min-disp, --span, --n-bins, --layer-key <element> (default layers/log1p).
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process SCANPY_HVG {
    tag id
    label 'process_gpu'
    container task.accelerator ?
        'quay.io/cellgeni/rapids-singlecell:0.18.1' :
        'community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90'

    input:
    record(id: String, adata: Path, log1p: Path?, counts_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, hvg: file("${output_prefix ?: id}.hvg.*"))

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
    scanpy_hvg.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${log1p ? "--log1p ${log1p}" : ''} \\
        ${counts_key ? "--counts-key ${counts_key}" : ''} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.hvg.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
