// INCLUDE MODULES
include { SEACellsAggregate } from './modules/SEAcells'
include { HierarchialAggregate } from './modules/hierarchial'
include { ToH5ad; AttachCellMetadata } from './modules/utils'

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
        --type                          Type of data to be aggregated (gex or atac)
        --celltype_label                Cell type label for the aggregation (required for hierarchial clustering; optional for SEACells)

        Optional options:
        --help                          Show this help message
        --output_dir                        Output directory for the aggregated files
        --raw                           Convert raw files to H5AD format (if raw files' paths are provided)
        --delimiter                     Specify delimiter if you want to add sample id to the barcode names
        --cell_metadata                 Specify metadata.csv file to attach to the .obs section of AnnData object
        --barcode_column                Barcode column name in the cell metadata file (specify "obs_names" if you want to use the .obs index column)

        SEACells options:
        --seacells.enabled              Enable SEACells aggregation
        --seacells.n_cells              Number of metacells to be aggregated
        --seacells.gamma                Gamma value for the adaptive bandwidth kernel
        --seacells.n_top_genes          Number of top variable genes to be used for aggregation (GEX only). Default: 2000
        --seacells.n_components         Number of components to be used for aggregation (PCA for GEX and LSI for ATAC). Default: 50
        --seacells.convergence_epsilon  Convergence epsilon for the optimization. Default: 0.00001
        --seacells.min_iterations       Minimum number of iterations for the optimization. Default: 10
        --seacells.max_iterations       Maximum number of iterations for the optimization. Default: 100
        --seacells.use_sparse           Use sparse matrix for the optimization. Default: false
        --seacells.precomputed          Use precomputed distance matrix for the aggregation.

        Hierarchial options:
        --hierarchial.enabled           Enable Hierarchial clustering aggregation
        --hierarchial.n_min             Minimum number of clusters for the aggregation
        --hierarchial.n_max             Maximum number of clusters for the aggregation
        --hierarchial.method            Method for the aggregation (kmeans or louvain)
        --hierarchial.n_top_genes       Number of top variable genes to be used for aggregation (GEX only). Default: 2000
        --hierarchial.n_components      Number of components to be used for aggregation (PCA for GEX and LSI for ATAC). Default: 50
        --hierarchial.n_neighbors       Number of neighbors for the aggregation. Default: 15
        --hierarchial.precomputed       Use precomputed distance matrix for the aggregation.
            
    Examples:
        1. Perform metacell aggregation using SEACells
            nextflow run main.nf --filelist filelist.csv --type gex --seacells.enabled --seacells.gamma 75

    == filelist.csv format ==
    item,filepath
    pbmc_10k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc_10k/pbmc10k/pbmc10k.h5ad
    pbmc_3k,/lustre/scratch126/cellgen/team361/data/pbmc/results/pbmc3k/pbmc3k/pbmc3k.h5ad
    ========================
    """.stripIndent()
}

workflow  {
    // Define the default parameters
    delimiter = params.delimiter
    barcode_column = params.barcode_column

    // Send help message if -help is used
    if ( params.help || !params.filelist ) {
        helpMessage()
        System.exit(0)
    }

    // Check if the filelist is provided
    if (!params.filelist || (!params.seacells.enabled && !params.hierarchial.enabled) || !params.type ) {
        log.error "Missing required arguments: filelist, type, seacells.enabled or hierarchial.enabled options"
        helpMessage()
        System.exit(1)
    }

    // Read the files from the filelist
    files = Channel.fromPath(params.filelist, checkIfExists: true)
                    .splitCsv(header: true, sep: ',')
                    .map { row -> tuple(row.item, row.filepath) }
    

    // Convert raw files to H5AD if needed
    if (params.raw) {
        // Convert H5 files to H5AD
        files = ToH5ad(
            files,
            delimiter ? delimiter : ""
        )

        // Set delimiter to empty string and barcode_column to "barcode" as it is done at for .h5ad file
        delimiter = ""
        barcode_column = "barcode"
    }

    // Attach cell metadata to .obs section of AnnData object if celltype annotation is provided
    if (params.cell_metadata) {

        // Check that barcode column is provided
        if (!barcode_column) {
            log.error "Cell metadata was provided, but --barcode_column is missing"
            System.exit(1)
        }
        files = AttachCellMetadata(
            files,
            params.cell_metadata,
            barcode_column
        )
    }
    

    // Run SEACells if enabled
    if (params.seacells.enabled) {
        // Check that all necessary arguments are provided
        if ( !params.seacells.n_cells && !params.seacells.gamma ) {
            log.error "Missing arguments for SEACells"
            helpMessage()
            System.exit(1)
        // Check that only one of n_cells or gamma is provided
        } else if (params.seacells.n_cells && params.seacells.gamma) {
            log.error "Both n_cells and gamma cannot be provided at the same time"
            helpMessage()
            System.exit(1)
        } else {
            // Run SEACells
            SEACellsAggregate(
                files,
                params.seacells.n_cells ? params.seacells.n_cells : "",
                params.seacells.gamma ? params.seacells.gamma : "",
                params.type,
                params.seacells.n_top_genes,
                params.seacells.n_components,
                params.celltype_label ? params.celltype_label : "",
                params.seacells.convergence_epsilon,
                params.seacells.min_iterations,
                params.seacells.max_iterations,
                params.seacells.use_sparse,
                params.seacells.precomputed ? params.hierarchial.precomputed : "",
                delimiter ? delimiter : ""
            )
        }
    }

    // Run Hierarchial clustering if enabled
    if (params.hierarchial.enabled) {
        // Check that all necessary arguments are provided
        if ( !params.hierarchial.n_min || !params.hierarchial.n_max || !params.celltype_label ) {
            log.error "Missing arguments for hierarchial clustering"
            helpMessage()
            System.exit(1)
        } else {
            // Run Hierarchial clustering
            HierarchialAggregate(
                files,
                params.hierarchial.n_min,
                params.hierarchial.n_max,
                params.type,
                params.hierarchial.n_top_genes,
                params.hierarchial.n_components,
                params.celltype_label,
                params.hierarchial.n_neighbors,
                params.hierarchial.precomputed ? params.hierarchial.precomputed : "",
                params.hierarchial.method,
                delimiter ? delimiter : ""
            )
        }
    }

}