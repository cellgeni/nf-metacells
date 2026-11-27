nextflow.enable.types = true

include { INGEST } from '../subworkflows/local/ingest'
include { REPRESENT_RNA } from '../subworkflows/local/represent_rna'
include { REPRESENT_ATAC } from '../subworkflows/local/represent_atac'
include { ANNDATA_ASSEMBLE } from '../modules/cellgeni/anndata/assemble'
include { SEACELLS } from '../modules/cellgeni/seacells'
include { HMETACELLS } from '../modules/cellgeni/hmetacells'
include { Sample; AssayState; StepOutput; Embedding; Membership; MethodFiles } from './types'

/*
 * Ingest every assay of each sample, compute the embedding of the sample's membership assay (the
 * --membership_assay when the sample has it, otherwise its only assay) unless --embedding_key names one in the
 * object, then run each selected method on that embedding: PCA for RNA, LSI or spectral (--atac_embedding) for ATAC.
 *
 * Until phase 5 (docs/architecture.md), SEACells and hmetacells are the v1 modules: they read the embedding from
 * obsm of the object they aggregate (--precomputed), so the step output is assembled with the base object first.
 * They see only the membership assay and write their own table formats.
 */
workflow METACELLS {
    take:
    samples: Channel<Sample>
    methods: List<String>
    membership_assay: String
    cell_metadata: Path?
    barcode_delimiter: String?
    embedding_key: String?
    atac_embedding: String
    rna_counts_key: String
    atac_counts_key: String

    main:
    ingest = INGEST(samples, cell_metadata, barcode_delimiter)

    bases = samples.map { s -> record(sample: s.id, basis: basisAssay(s, membership_assay)) }
    basis = ingest.states.join(bases, by: 'sample').filter { st -> st.assay == st.basis }

    compute = basis.filter { _st -> !embedding_key }
    rna = REPRESENT_RNA(compute.filter { st -> st.assay == 'rna' }, rna_counts_key)
    atac = REPRESENT_ATAC(compute.filter { st -> st.assay == 'atac' }, atac_embedding, atac_counts_key)
    computed = compute.join(
        rna.embeddings.mix(atac.embeddings)
            .map { o -> record(id: "${o.sample}.${o.assay}", step: o.step, output: o.output) },
        by: 'id'
    )
    assembled = ANNDATA_ASSEMBLE(
        computed.map { st -> record(id: st.id, adata: st.adata, sources: [st.output], output_prefix: null) }
    )
    embeddings = computed.join(assembled, by: 'id')
        .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: "X_${st.step}", adata: st.assembled) }
        .mix(
            basis.filter { _st -> embedding_key }
                .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: embedding_key, adata: st.adata) }
        )
    v1_input = embeddings.map { e -> tuple([id: e.sample, assay: e.assay, embedding_key: e.key], e.adata) }

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
    representations: Channel<StepOutput> = rna.steps.mix(atac.steps)
    embeddings: Channel<Embedding> = embeddings
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
