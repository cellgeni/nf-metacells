nextflow.enable.types = true

include { ADATA_CONVERT } from '../../../modules/cellgeni/adata/convert'
include { ADATA_ATTACHOBS } from '../../../modules/cellgeni/adata/attachobs'
include { Sample; AssayState } from '../../../workflows/types'

/*
 * One AnnData per assay of each sample. 10x inputs (an .h5 file or an mtx folder) are converted with adata-cli;
 * AnnData inputs (.h5ad, .zarr) are used as given. Cell metadata is then attached and/or barcodes are made unique
 * per sample when asked; a sample's own metadata file wins over the global one.
 *
 * `written` holds the states whose AnnData this subworkflow wrote (converted or attached), for publishing.
 */
workflow INGEST {
    take:
    samples: Channel<Sample>
    cell_metadata: Path?
    barcode_delimiter: String?

    main:
    rna = samples.filter { s -> s.rna }
        .map { s -> record(id: "${s.id}.rna", sample: s.id, assay: 'rna', source: s.rna, metadata: s.metadata ?: cell_metadata) }
    atac = samples.filter { s -> s.atac }
        .map { s -> record(id: "${s.id}.atac", sample: s.id, assay: 'atac', source: s.atac, metadata: s.metadata ?: cell_metadata) }
    assays = rna.mix(atac)

    // 10x inputs → AnnData
    tenx = assays.filter { a -> !isAnnData(a.source) }
    converted = ADATA_CONVERT(
        tenx.map { a -> record(id: a.id, source: a.source, convert_args: convertArgs(a.assay, a.source), output_prefix: a.id, meta: null) }
    )
    stored = assays.filter { a -> isAnnData(a.source) }
        .map { a -> record(id: a.id, sample: a.sample, assay: a.assay, adata: a.source, metadata: a.metadata, written: false) }
        .mix(
            tenx.join(converted, by: 'id')
                .map { a -> record(id: a.id, sample: a.sample, assay: a.assay, adata: a.adata, metadata: a.metadata, written: true) }
        )

    // Cell metadata and/or barcode suffix
    attach = stored.filter { st -> st.metadata || barcode_delimiter }
    attached = ADATA_ATTACHOBS(
        attach.map { st -> record(id: st.id, adata: st.adata, metadata: st.metadata, sample_id: st.sample, output_prefix: st.id, meta: null) }
    )
    ready = stored.filter { st -> !(st.metadata || barcode_delimiter) }
        .mix(
            attach.join(attached, by: 'id')
                .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, adata: st.adata, metadata: st.metadata, written: true) }
        )

    states = ready.map { st -> record(id: st.id, sample: st.sample, assay: st.assay, adata: st.adata) }
    written = ready.filter { st -> st.written }
        .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, adata: st.adata) }

    emit:
    states: Channel<AssayState> = states
    written: Channel<AssayState> = written
}

def isAnnData(path: Path) -> Boolean {
    return path.name.endsWith('.h5ad') || path.name.endsWith('.zarr')
}

// --feature-type only for inputs that record feature types, guessed from the path until cellgeni/adata-cli#127 is
// fixed (a Cell Ranger 2 .h5 still fails). RNA var names are gene symbols, made unique, as scanpy reads 10x data.
def convertArgs(assay: String, source: Path) -> String {
    def typed = source.isDirectory()
        ? source.resolve('features.tsv.gz').exists() || source.resolve('features.tsv').exists()
        : source.name.endsWith('.h5')
    def featureType = assay == 'atac' ? 'Peaks' : 'Gene Expression'
    return (typed ? "--feature-type '${featureType}' " : '') + (assay == 'rna' ? '--var-names gene_symbols --make-unique' : '')
}
