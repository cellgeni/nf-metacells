// INCLUDE MODULES
include { SEACELLS } from './modules/cellgeni/seacells'
include { HMETACELLS } from './modules/cellgeni/hmetacells'
include { ADATA_CONVERT } from './modules/cellgeni/adata/convert'
include { ADATA_ATTACHOBS } from './modules/cellgeni/adata/attachobs'

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
        --raw                           Convert 10x .h5 files or mtx folders to H5AD with adata-cli
        --delimiter                     Rename barcodes to <barcode><delimiter><id> (originals kept in obs['barcode'])
        --cell_metadata                 CSV (barcode, sample_id, ...) to attach to .obs; cells without metadata are dropped
        --barcode_column                obs column matched to the metadata barcodes. Default: obs_names
        --allow_dtype_change            Accept digit-only categorical obs columns becoming integers when attaching metadata

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
    files = channel.fromPath(params.filelist, checkIfExists: true)
                    .splitCsv(header: true, sep: ',')
                    .map { row -> tuple([id: row.id], file(row.path, checkIfExists: true)) }

    // Convert raw 10x files to H5AD with adata-cli; --feature-type only for inputs that record feature types,
    // guessed from the path until cellgeni/adata-cli#127 is fixed (a Cell Ranger 2 .h5 still fails)
    if (params.raw) {
        def feature_type = params.type == 'atac' ? 'Peaks' : 'Gene Expression'
        converted = ADATA_CONVERT(
            files.map { meta, path ->
                def typed = path.isDirectory()
                    ? (path.resolve('features.tsv.gz').exists() || path.resolve('features.tsv').exists())
                    : path.name.endsWith('.h5')
                def args = (typed ? "--feature-type '${feature_type}' " : '') + (params.type == 'gex' ? '--var-names gene_symbols --make-unique' : '')
                record(id: meta.id, source: path, convert_args: args, output_prefix: null, meta: null)
            }
        )
        files = converted.map { r -> tuple([id: r.id], r.adata) }
    }

    // Attach cell metadata and/or make barcodes unique per sample
    if (params.cell_metadata || params.delimiter) {
        def metadata = params.cell_metadata ? file(params.cell_metadata, checkIfExists: true) : null
        attached = ADATA_ATTACHOBS(
            files.map { meta, path ->
                record(id: meta.id, adata: path, metadata: metadata, sample_id: meta.id, output_prefix: null, meta: null)
            }
        )
        files = attached.map { r -> tuple([id: r.id], r.adata) }
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
            SEACELLS(
                files
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
            HMETACELLS(files)
        }
    }

}
