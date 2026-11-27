nextflow.enable.types = true

/*
 * Latent semantic indexing of a TF-IDF layer into obsm/X_lsi, varm/LSI and uns/lsi, as muon's lsi does: SciPy
 * svds on CPU, cupyx svds when the task has an accelerator. ext.args: --n-comps <int> (default 50), --matrix-key
 * (default layers/tfidf), --no-scale-embeddings, --drop-first. Components are signed with sklearn's svd_flip
 * convention.
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process ATAC_LSI {
    tag id
    label 'process_gpu'
    container task.accelerator ? 'quay.io/cellgeni/rapids-singlecell:0.18.1' : 'quay.io/cellgeni/scanpy:1.12.6'

    input:
    record(id: String, adata: Path, tfidf: Path?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, lsi: file("${output_prefix ?: id}.lsi.*"))

    topic:
    record(
        process: task.process,
        tool: task.accelerator ? 'cupy' : 'scipy',
        version: eval("scio.py version ${task.accelerator ? 'cupy-cuda12x' : 'scipy'}")
    ) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    atac_lsi.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${tfidf ? "--tfidf ${tfidf}" : ''} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.lsi.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
