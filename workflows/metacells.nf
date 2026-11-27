nextflow.enable.types = true

include { INGEST } from '../subworkflows/local/ingest'
include { SEACELLS } from '../modules/cellgeni/seacells'
include { HMETACELLS } from '../modules/cellgeni/hmetacells'
include { Sample; AssayState; Membership; MethodFiles } from './types'

/*
 * Ingest every assay of each sample, then run each selected method on the sample's membership assay: the
 * --membership_assay when the sample has it, otherwise its only assay.
 *
 * Until the per-step modules land (docs/architecture.md, phases 4-5), SEACells and hmetacells are the v1 modules:
 * they preprocess internally, see only the membership assay, and write their own table formats.
 */
workflow METACELLS {
    take:
    samples: Channel<Sample>
    methods: List<String>
    membership_assay: String
    cell_metadata: Path?
    barcode_delimiter: String?

    main:
    ingest = INGEST(samples, cell_metadata, barcode_delimiter)

    bases = samples.map { s -> record(sample: s.id, basis: basisAssay(s, membership_assay)) }
    basis = ingest.states.join(bases, by: 'sample').filter { st -> st.assay == st.basis }
    v1_input = basis.map { st -> tuple([id: st.sample, assay: st.assay], st.adata) }

    seacells = SEACELLS(v1_input.filter { _st -> 'seacells' in methods })
    hmetacells = HMETACELLS(v1_input.filter { _st -> 'hmetacells' in methods })

    memberships = seacells.csv
        .map { meta, tables -> membership(meta, 'seacells', asList(tables).find { f -> isSeacellsHard(f) } ?: asList(tables)[0]) }
        .mix(hmetacells.csv.map { meta, table -> membership(meta, 'hmetacells', table) })

    extras = seacells.h5ad.mix(seacells.csv).mix(seacells.npy).mix(seacells.pdf).mix(seacells.pkl)
        .map { meta, files -> methodFiles(meta, 'seacells', asList(files).findAll { f -> !isSeacellsHard(f) }) }
        .filter { r -> r.files }

    emit:
    states: Channel<AssayState> = ingest.states
    written: Channel<AssayState> = ingest.written
    memberships: Channel<Membership> = memberships
    extras: Channel<MethodFiles> = extras
}

// The v1 SEACells hard assignments serve as its membership table (a stub run writes one table of another name)
def isSeacellsHard(file: Path) -> Boolean {
    return file.name == 'seacell_metacell_assignments.csv'
}

def basisAssay(s: Sample, preferred: String) -> String {
    if (preferred == 'atac') {
        return s.atac ? 'atac' : 'rna'
    }
    return s.rna ? 'rna' : 'atac'
}

def membership(meta: Map, method: String, table: Path) -> Membership {
    return record(id: "${meta.id}.${meta.assay}.${method}", sample: meta.id, assay: meta.assay, method: method, membership: table)
}

def methodFiles(meta: Map, method: String, files: List<Path>) -> MethodFiles {
    return record(id: "${meta.id}.${meta.assay}.${method}.${files ? files[0].name : 'none'}", sample: meta.id, method: method, files: files)
}

// A legacy path output is a single path when its glob matches one file and a list otherwise
def asList(files) -> List<Path> {
    return files instanceof List ? files : [files]
}
