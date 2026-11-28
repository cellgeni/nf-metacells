nextflow.enable.types = true

/*
 * SEACells metacells from a precomputed embedding: `<prefix>.membership.csv` (barcode,metacell) and
 * `<prefix>.seacells.json` (convergence, RSS, metacell sizes, versions). The unified SEACells model runs on CPU
 * and, when the task has an accelerator, on the device; ext.args `--use-sparse` asks for the legacy sparse CPU
 * backend instead. One of `--n-metacells <int>` or `--gamma <int>` is required (ext.args); the obsm key comes
 * from the record. Other ext.args: --n-waypoint-eigs, --n-neighbors, --convergence-epsilon, --min-iter,
 * --max-iter, --use-sparse, --seed.
 * Non-convergence is reported in the JSON rather than raised: the labels are set before SEACells complains, and
 * with a seed a retry reproduces the same partition.
 */
process SEACELLS {
    tag id
    label 'process_gpu'
    container task.accelerator ?
        'quay.io/cellgeni/seacells:a1173c8-gpu' :
        'quay.io/cellgeni/seacells:a1173c8-cpu'

    input:
    record(id: String, adata: Path, embedding: Path?, embedding_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(
        id: id,
        membership: file("${output_prefix ?: id}.membership.csv"),
        summary: file("${output_prefix ?: id}.seacells.json")
    )

    topic:
    record(process: task.process, tool: 'SEACells', version: eval('scio.py version SEACells')) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    export PYTHONNOUSERSITE=1 NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}/numba MPLCONFIGDIR=\${TMPDIR:-/tmp}/matplotlib
    seacells_membership.py \\
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
    printf '{"method": "seacells", "converged": true}\\n' > ${output_prefix ?: id}.seacells.json
    """
}
