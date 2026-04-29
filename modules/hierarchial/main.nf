process HierarchialAggregate {
    tag "Running hierarchial clustering for ${meta.id}"
    container 'quay.io/cellgeni/metacells-python:latest'
    
    input:
    tuple val(meta), path(adata)
    val n_min
    val n_max
    val celltype_label

    output:
    path("hierarchial_metacells.csv")
    
    script:
    def args = task.ext.args ?: "--n_min 5 --n_max 20 --method ward --n_top_genes 2000 --n_components 50 --n_neighbors 15"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hierarchial_metacells.py \
        ${args} \
        --adata ${adata} \
        --sample ${prefix} \
        --n_min ${n_min} \
        --n_max ${n_max} \
        --celltype_label ${celltype_label} \
        --output hierarchial_metacells.csv 
    """
}