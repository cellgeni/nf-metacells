nextflow.enable.types = true

/*
 * nf-metacells: metacells (SEACells, hierarchical) from scRNA-seq and scATAC-seq AnnData.
 * Design, contracts and roadmap: docs/architecture.md
 */

include { paramsHelp; paramsSummaryLog; samplesheetToList; validateParameters } from 'plugin/nf-schema'
include { METACELLS } from './workflows/metacells'
include { AssayState; Membership; MethodFiles } from './workflows/types'

// Documented in nextflow_schema.json; `outdir` and `publish_mode` are in nextflow.config, which uses them
params {
    // Input
    input: Path?
    cell_metadata: Path?
    barcode_column: String = 'obs_names'
    barcode_delimiter: String?
    allow_dtype_change: Boolean = false

    // Methods
    methods: String = 'seacells'
    membership_assay: String = 'rna'
    metacell_size: Integer = 75
    constrain_by: String?
    eval_labels: String?

    // Representation (computed inside the methods until the per-step modules land)
    embedding_key: String?
    hvg_n_top_genes: Integer = 2000
    pca_n_comps: Integer = 50
    lsi_n_comps: Integer = 50

    // SEACells
    seacells_n_metacells: Integer?
    seacells_n_waypoint_eigs: Integer = 10
    seacells_convergence_epsilon: Float = 0.00001
    seacells_min_iter: Integer = 10
    seacells_max_iter: Integer = 50
    seacells_use_sparse: Boolean = false

    // Hierarchical metacells
    hmetacells_n_min: Integer?
    hmetacells_n_max: Integer?
    hmetacells_method: String = 'louvain'
    hmetacells_n_neighbors: Integer = 15

    // GPU (read by conf/base.config for processes labelled process_gpu)
    gpu_steps: String = 'none'
    gpuqueue: String = 'gpu'
    gpu_gmem: Integer = 40000
    costcode: String?

    // --help, --help <param>, --helpFull
    help: String?
    helpFull: Boolean = false
}

workflow {
    main:
    if (params.help || params.helpFull) {
        exit(0, paramsHelp(
            parameter: !params.help || params.help == 'true' ? null : params.help,
            beforeText: helpBefore(),
            command: 'nextflow run cellgeni/nf-metacells --input samplesheet.csv --methods seacells'
        ))
    }
    validateParameters()
    log.info(paramsSummaryLog())

    def methods = params.methods.tokenize(',')
    if ('hmetacells' in methods && !(params.hmetacells_n_min && params.hmetacells_n_max && params.constrain_by)) {
        error('hmetacells needs --hmetacells_n_min, --hmetacells_n_max and --constrain_by (the obs column to build metacells within)')
    }

    // Rows come in the property order of assets/schema_input.json; empty cells come as []
    samples = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { row -> record(id: row[0], rna: row[1] ?: null, atac: row[2] ?: null, metadata: row[3] ?: null) }

    run = METACELLS(samples, methods, params.membership_assay, params.cell_metadata, params.barcode_delimiter)

    publish:
    ingested = run.written
    memberships = run.memberships
    method_files = run.extras
    versions = channel.topic('versions').unique()
}

output {
    ingested: Channel<AssayState> {
        path { st -> "${st.sample}/ingest" }
        index {
            path 'ingested.csv'
            header true
        }
    }
    memberships: Channel<Membership> {
        path { m -> "${m.sample}/${m.method}" }
        index {
            path 'memberships.csv'
            header true
        }
    }
    method_files: Channel<MethodFiles> {
        path { f -> "${f.sample}/${f.method}" }
    }
    versions {
        index {
            path 'pipeline_info/versions.csv'
            header true
        }
    }
}

def helpBefore() -> String {
    return '''
    =============
    nf-metacells
    =============
    Metacells from scRNA-seq and scATAC-seq AnnData: SEACells and hierarchical (per-group) metacells.

    == samplesheet (comma-separated) ==
    id,rna,atac,metadata
    S1,/path/S1_gex.h5ad,/path/S1_atac.h5ad,
    S2,/path/S2/filtered_feature_bc_matrix.h5,/path/S2/filtered_feature_bc_matrix.h5,/path/S2_metadata.csv
    S3,/path/S3/filtered_feature_bc_matrix,,

    rna and atac take an .h5ad or .zarr store, a 10x .h5 file or a 10x mtx folder (a multiome .h5 may be given for
    both); every sample needs at least one of them. metadata is an optional CSV (barcode, annotations...) for that
    sample; --cell_metadata gives one for all samples (barcode, sample_id, annotations...).
    '''.stripIndent()
}
