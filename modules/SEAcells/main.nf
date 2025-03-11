// The module that contains the processes to run SEAcells metacells aggregation
process SEACellsAggregate {
    tag "Running SEAcells metacells aggregation foor ${item}"
    input:
        tuple val(item), path(adata)
        val(n_cells)
        val(gamma)
        val(type)
        val(n_top_genes)
        val(n_components)
        val(celltype_label)
        val(convergence_epsilon)
        val(min_iterations)
        val(max_iterations)
        val(use_sparse)
    output:
        path("${item}")
    script:
        """
        mkdir -p ${item}
        seacells_aggregate.py \
            --adata ${adata} \
            ${n_cells ? "--n_metacells ${n_cells}" : ""} \
            ${gamma ? "--gamma ${gamma}" : ""} \
            --type ${type} \
            --output_dir ${item} \
            --n_top_genes ${n_top_genes} \
            --n_components ${n_components} \
            ${celltype_label ? "--celltype_label ${celltype_label}" : ""} \
            --convergence_epsilon ${convergence_epsilon} \
            --min_iter ${min_iterations} \
            --max_iter ${max_iterations} \
            ${use_sparse ? "--use_sparse" : ""}
        """
}