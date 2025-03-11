// IMPORT MODULES
include { SEACellsAggregate } from './modules/SEAcells/main.nf'

// HELP MESSAGE
def helpMessage() {
    log.info """
    Usage: nextflow run main.nf --input_dir <input_dir> --output_dir <output_dir> --n_cells <n_cells> --gamma <gamma> --type <type> --n_top_genes <n_top_genes> --n_components <n_components> --celltype_label <celltype_label> --convergence_epsilon <convergence_epsilon> --min_iterations <min_iterations> --max_iterations <max_iterations>
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