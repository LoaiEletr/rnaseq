//
// Prepare reference genome and build alignment indices for RNA-seq analysis
// Processes: Genome decompression, GTF annotation processing, alignment index construction
// Builds indices for: HISAT2 (genome alignment), Kallisto/Salmon (pseudo-alignment),
// SortMeRNA (rRNA filtering), BBMap bbsplit (contaminant removal), GXF2BED for RSeQC
// Conditional execution based on workflow configuration parameters
//

include {
    GUNZIP as GUNZIP_GENOME ;
    GUNZIP as GUNZIP_GTF
} from '../../../modules/local/gunzip'
include { UNTAR } from '../../../modules/local/untar'
include { GXF2BED } from '../../../modules/local/gxf2bed'
include { KALLISTO_INDEX } from '../../../modules/local/kallisto/index'
include { SALMON_INDEX } from '../../../modules/local/salmon/index'
include { HISAT2_EXTRACTEXONS } from '../../../modules/local/hisat2/extractexons'
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/local/hisat2/extractsplicesites'
include { HISAT2_BUILD } from '../../../modules/local/hisat2/build'
include { SORTMERNA as SORTMERNA_INDEX } from '../../../modules/local/sortmerna'
include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_INDEX } from '../../../modules/local/bbmap/bbsplit'

workflow PREPARE_GENOME {
    take:
    val_fasta // string: Path to reference genome FASTA (.gz)
    val_transcriptome // string: Path to transcriptome FASTA (.gz) for pseudo-alignment
    val_gtf // string: Path to gene annotation GTF (.gz)
    val_gtf_isoform // string: Path to extended GTF (.gz) for isoform analysis
    val_contaminant_fasta // string: Path to contaminant genome FASTA (.gz) for BBSplit
    val_rrna_db // string: Path to rRNA database file (.tar.gz) for SortMeRNA
    val_rrna_db_type // string: default, fast, or sensitive

    main:
    ch_versions = channel.empty()
    ch_fasta_uncompressed = channel.empty()
    ch_gtf_uncompressed = channel.empty()
    ch_rrna_db_fasta = channel.empty()
    ch_hisat2_index = channel.empty()
    ch_kallisto_index = channel.empty()
    ch_salmon_index = channel.empty()
    ch_sortmerna_index = channel.empty()
    ch_bbsplit_index = channel.empty()
    ch_bed = channel.empty()

    // Convert string paths to file channels
    ch_fasta = params.aligner == "hisat2" || !params.skip_bbsplit ? channel.fromPath(val_fasta, checkIfExists: true) : channel.empty()
    ch_transcriptome = val_transcriptome ? channel.fromPath(val_transcriptome, checkIfExists: true) : channel.empty()
    ch_gtf = channel.fromPath(val_gtf, checkIfExists: true)
    ch_gtf_isoform = "DIU" in params.analysis_method.split(",") ? channel.fromPath(val_gtf_isoform, checkIfExists: true) : channel.empty()
    ch_contaminant_fasta = val_contaminant_fasta && !params.skip_bbsplit ? channel.fromPath(val_contaminant_fasta, checkIfExists: true) : channel.empty()
    ch_rrna_db = !params.skip_sortmerna ? channel.fromPath(val_rrna_db, checkIfExists: true) : channel.empty()

    // 1. HISAT2 genome alignment index construction (if aligner is hisat2)
    if (params.aligner == "hisat2") {
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

        // Extract splice sites from GTF
        HISAT2_EXTRACTSPLICESITES(ch_gtf_uncompressed)
        ch_splice_sites = HISAT2_EXTRACTSPLICESITES.out.splice_sites
        ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

        // Extract exons from GTF
        HISAT2_EXTRACTEXONS(ch_gtf_uncompressed)
        ch_exon_sites = HISAT2_EXTRACTEXONS.out.exons_sites
        ch_versions = ch_versions.mix(HISAT2_EXTRACTEXONS.out.versions)

        // Build HISAT2 index with reference genome, splice sites, and exons
        HISAT2_BUILD(
            ch_fasta_uncompressed,
            ch_splice_sites,
            ch_exon_sites,
        )
        ch_hisat2_index = HISAT2_BUILD.out.index
        ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)

        // Convert GTF to BED format for RSeQC (if RSeQC modules are enabled)
        if (params.rseqc_modules) {
            GXF2BED(
                ch_gtf
            )
            ch_bed = GXF2BED.out.bed
            ch_versions = ch_versions.mix(GXF2BED.out.versions)
        }
    }

    // 2. Kallisto transcriptome index construction (if pseudo_aligner is kallisto)
    if (params.pseudo_aligner == "kallisto") {
        KALLISTO_INDEX(
            ch_transcriptome
        )
        ch_kallisto_index = KALLISTO_INDEX.out.index
        ch_versions = ch_versions.mix(KALLISTO_INDEX.out.versions)
    }

    // 3. Salmon transcriptome index construction
    SALMON_INDEX(
        ch_transcriptome
    )
    ch_salmon_index = SALMON_INDEX.out.index
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    // 4. SortMeRNA rRNA database index construction (if not skipped)
    if (!params.skip_sortmerna) {
        UNTAR(ch_rrna_db)
        ch_rrna_db_fasta = UNTAR.out.untar
        ch_versions = ch_versions.mix(UNTAR.out.versions)

        SORTMERNA_INDEX(
            [[], []],
            [],
            ch_rrna_db_fasta,
            val_rrna_db_type,
            true,
        )
        ch_sortmerna_index = SORTMERNA_INDEX.out.index
        ch_versions = ch_versions.mix(SORTMERNA_INDEX.out.versions)
    }

    // 5. BBMap bbsplit contamination index construction (if not skipped)
    if (!params.skip_bbsplit) {
        BBMAP_BBSPLIT_INDEX(
            [[id: 'index'], []],
            [],
            ch_fasta,
            ch_contaminant_fasta,
            true,
        )
        ch_bbsplit_index = BBMAP_BBSPLIT_INDEX.out.index
        ch_versions = ch_versions.mix(BBMAP_BBSPLIT_INDEX.out.versions)
    }

    emit:
    fasta_compressed = ch_fasta // channel: [ genome.fasta.gz ]
    fasta_uncompressed = ch_fasta_uncompressed // channel: [ genome.fasta ]
    transcriptome = ch_transcriptome // channel: [ transcriptome.fasta.gz ]
    gtf_compressed = ch_gtf // channel: [ genes.gtf.gz ]
    gtf_uncompressed = ch_gtf_uncompressed // channel: [ genes.gtf ]
    gtf_isoform = ch_gtf_isoform // channel: [ genes_with_patches.gtf.gz ]
    contaminant_fasta = ch_contaminant_fasta // channel: [ contaminant.fasta.gz ]
    rrna_db_fasta = ch_rrna_db_fasta // channel: [ rrna.fasta ]
    bed = ch_bed // channel: [ genes.bed ]
    hisat2_index = ch_hisat2_index // channel: [ hisat2_index_files ]
    kallisto_index = ch_kallisto_index // channel: [ kallisto_index_files ]
    salmon_index = ch_salmon_index // channel: [ salmon_index_files ]
    sortmerna_index = ch_sortmerna_index // channel: [ sortmerna_index_files ]
    bbsplit_index = ch_bbsplit_index // channel: [ bbsplit_index_files ]
    versions = ch_versions // channel: [ versions.yml ]
}
