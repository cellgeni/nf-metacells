// The module that contains the processes to run SEAcells metacells aggregation
process SEACellsAggregate {
    tag "Running SEAcells metacells aggregation for ${meta.id}"
    container 'docker://quay.io/cellgeni/seacells:latest'

    input:
    tuple val(meta), path(adata)

    output:
    path("*")

    script:
    def args = task.ext.args ?: "--type gex --n_top_genes 2000 --n_components 50 --convergence_epsilon 0.00001 --min_iter 10 --max_iter 50"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seacells_aggregate.py \
        ${args} \
        --adata ${adata} \
        --sample ${prefix} \
        --output_dir .
    """
}