// IMPORT MODULES
include { SEACellsAggregate } from './modules/SEAcells/main.nf'

// HELP MESSAGE
def helpMessage() {
    log.info """
    ===========================
    Pipeline to aggregate metacells/pseudobulks
    ===========================
    This pipeline performs metacell-aggregation (only SEACells option is available at the moment)
    Usage: nextflow run main.nf [OPTIONS]
        Required options:
        --filelist                      File containing the list of files to be aggregated

        SEACells options:
        --seacells.enabled              Enable SEACells aggregation
        --seacells.n_cells              Number of metacells to be aggregated
        --seacells.gamma                Gamma value for the adaptive bandwidth kernel
        --seacells.type                 Type of aggregation (gex or atac)
        --seacells.n_top_genes          Number of top variable genes to be used for aggregation (GEX only)
        --seacells.n_components         Number of components to be used for aggregation (PCA for GEX and LSI for ATAC)
        --seacells.celltype_label       Label for the celltype column in the metadata (used to calculate aggregation metrics)
        --seacells.convergence_epsilon  Convergence epsilon for the optimization
        --seacells.min_iterations       Minimum number of iterations for the optimization
        --seacells.max_iterations       Maximum number of iterations for the optimization
        --seacells.use_sparse           Use sparse matrix for the optimization
            
    Examples:
        1. Perform metacell aggregation using SEACells
            nextflow run main.nf --filelist filelist.csv --seacells.enabled --seacells.gamma 75 --seacells.type gex

    == filelist.csv format ==
    item,filepath
    pbmc_10k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc_10k/pbmc10k/pbmc10k.h5ad
    pbmc_3k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc3k/pbmc3k/pbmc3k.h5ad
    ========================
    """.stripIndent()
}

workflow  {
    // Send help message if no arguments are provided or -help is used
    if (params.help || !params.filelist || !params.seacells.enabled) {
        helpMessage()
        System.exit(0)
    }

    // Run SEACells if enabled
    if (params.seacells.enabled) {
        // Check that all necessary arguments are provided
        if ( ( !params.seacells.n_cells && !params.seacells.gamma ) || !params.seacells.type ) {
            log.error "Missing arguments for SEACells"
            helpMessage()
            System.exit(1)
        } else if (params.seacells.n_cells && params.seacells.gamma) {
            log.error "Both n_cells and gamma cannot be provided at the same time"
            helpMessage()
            System.exit(1)
        } else {
            // Read the files from the filelist
            files = Channel.fromPath(params.filelist, checkIfExists: true)
                           .splitCsv(header: true, sep: ',')
                            .map { row -> tuple(row.item, row.filepath) }
            // Run SEACells
            SEACellsAggregate(
                files,
                params.seacells.n_cells ? params.seacells.n_cells : "",
                params.seacells.gamma ? params.seacells.gamma : "",
                params.seacells.type,
                params.seacells.n_top_genes,
                params.seacells.n_components,
                params.seacells.celltype_label ? params.seacells.celltype_label : "",
                params.seacells.convergence_epsilon,
                params.seacells.min_iterations,
                params.seacells.max_iterations,
                params.seacells.use_sparse
            )
        }
    }

}