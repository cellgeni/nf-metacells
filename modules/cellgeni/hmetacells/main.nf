nextflow.enable.types = true

/*
 * Hierarchical metacells built within the groups of an obs column: `<prefix>.membership.csv` (barcode,metacell)
 * and `<prefix>.hmetacells.json` (metacell sizes, groups, unassigned cells, versions). Each group's kNN graph of
 * the embedding is cut by a scikit-network hierarchy into balanced clusters of at most `--n-max` cells; the graph
 * is built with scanpy on CPU and with rapids-singlecell when the task has an accelerator. ext.args: `--constrain-by
 * <obs column>`, `--n-min <int>` and `--n-max <int>` are required; also --method louvain|paris, --n-neighbors
 * (default 15), --n-dims and --seed (default 0). The obsm key comes from the record.
 */
process HMETACELLS {
    tag id
    label 'process_gpu'
    container task.accelerator ?
        'quay.io/cellgeni/rapids-singlecell:0.18.1' :
        'community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90'

    input:
    record(id: String, adata: Path, embedding: Path?, embedding_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(
        id: id,
        membership: file("${output_prefix ?: id}.membership.csv"),
        summary: file("${output_prefix ?: id}.hmetacells.json")
    )

    topic:
    record(process: task.process, tool: 'scikit-network', version: eval('scio.py version scikit-network')) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    export PYTHONNOUSERSITE=1 NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}/numba MPLCONFIGDIR=\${TMPDIR:-/tmp}/matplotlib
    hmetacells_membership.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${embedding ? "--embedding ${embedding}" : ''} \\
        ${embedding_key ? "--embedding-key ${embedding_key}" : ''} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    """
    printf 'barcode,metacell\\n' > ${output_prefix ?: id}.membership.csv
    printf '{"method": "hmetacells", "n_metacells": 0}\\n' > ${output_prefix ?: id}.hmetacells.json
    """
}
