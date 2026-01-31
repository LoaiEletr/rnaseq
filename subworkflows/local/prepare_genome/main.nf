//
// Decompress and prepare reference genome files for RNA-seq analysis
// Processes: Genome FASTA decompression, GTF annotation decompression, rRNA database extraction, contaminant genome preparation
// Conditional handling: Extended GTF for isoform analysis (DIU only), rRNA database for SortMeRNA, contaminant genomes for BBSplit
//

include {
    GUNZIP as GUNZIP_GENOME ;
    GUNZIP as GUNZIP_GTF
} from '../../../modules/local/gunzip'
include { UNTAR } from '../../../modules/local/untar'

workflow PREPARE_GENOME {
    take:
    val_fasta // string: Path to reference genome FASTA (.gz)
    val_transcriptome // string: Path to transcriptome FASTA (.gz) for pseudo-alignment
    val_gtf // string: Path to gene annotation GTF (.gz)
    val_gtf_isoform // string: Path to extended GTF (.gz) for isoform analysis
    val_contaminant_fasta // string: Path to contaminant genome FASTA (.gz) for BBSplit
    val_rrna_db // string: Path to rRNA database file (.tar.gz) for SortMeRNA

    main:
    ch_versions = channel.empty()

    // Convert string paths to file channels
    ch_fasta = channel.fromPath(val_fasta, checkIfExists: true)
    ch_transcriptome = channel.fromPath(val_transcriptome, checkIfExists: true)
    ch_gtf = channel.fromPath(val_gtf, checkIfExists: true)
    ch_gtf_isoform = "DIU" in params.analysis_method.split(",") ? channel.fromPath(val_gtf_isoform, checkIfExists: true) : channel.empty()
    ch_contaminant_fasta = val_contaminant_fasta ? channel.fromPath(val_contaminant_fasta, checkIfExists: true) : channel.empty()
    ch_rrna_db = !params.skip_sortmerna ? channel.fromPath(val_rrna_db, checkIfExists: true) : channel.empty()

    // Uncompress genome FASTA
    GUNZIP_GENOME(ch_fasta.map { it -> [[id: it.baseName], it] }.collect())
    ch_fasta_uncompressed = GUNZIP_GENOME.out.gunzip.map { basename, fasta ->
        fasta
    }
    ch_versions = ch_versions.mix(GUNZIP_GENOME.out.versions)

    // Uncompress gene annotation GTF
    GUNZIP_GTF(ch_gtf.map { it -> [[id: it.baseName], it] }.collect())
    ch_gtf_uncompressed = GUNZIP_GTF.out.gunzip.map { basename, gtf ->
        gtf
    }
    ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

    ch_rrna_db_fasta = channel.empty()
    // Uncompress rRNA database
    if (!params.skip_sortmerna) {
        UNTAR(ch_rrna_db)
        ch_rrna_db_fasta = UNTAR.out.untar
        ch_versions = ch_versions.mix(UNTAR.out.versions)
    }

    emit:
    fasta_compressed = ch_fasta // channel: [ path(genome.fasta.gz) ]
    fasta_uncompressed = ch_fasta_uncompressed // channel: [ path(genome.fasta) ]
    transcriptome = ch_transcriptome // channel: [ path(transcriptome.fasta.gz) ]
    gtf_compressed = ch_gtf // channel: [ path(genes.gtf.gz) ]
    gtf_uncompressed = ch_gtf_uncompressed // channel: [ path(genes.gtf) ]
    gtf_isoform = ch_gtf_isoform // channel: [ path(genes_with_patches.gtf.gz) ]
    contaminant_fasta = ch_contaminant_fasta // channel: [ path(contaminant.fasta.gz) ]
    rrna_db_fasta = ch_rrna_db_fasta // channel: [ path(rrna.fasta) ]
    versions = ch_versions // channel: [ versions.yml ]
}
