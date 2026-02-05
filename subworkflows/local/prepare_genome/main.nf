//
// Prepare reference genome and build alignment indices for RNA-seq analysis
// Processes: Genome decompression, GTF annotation processing, alignment index construction
// Builds indices for: HISAT2 (genome alignment), Kallisto/Salmon (pseudo-alignment),
// SortMeRNA (rRNA filtering), BBMap bbsplit (contaminant removal), GXF2BED for RSeQC
// Conditional execution based on workflow configuration parameters
//

include {
    GUNZIP as GUNZIP_GENOME ;
    GUNZIP as GUNZIP_GTF ;
    GUNZIP as GUNZIP_KALLISTO_INDEX
} from '../../../modules/local/gunzip'
include {
    UNTAR as UNTAR_BBSPLIT_INDEX ;
    UNTAR as UNTAR_SORTMERNA_INDEX ;
    UNTAR as UNTAR_HISAT2_INDEX ;
    UNTAR as UNTAR_SALMON_INDEX
} from '../../../modules/local/untar'

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
    val_fasta // string: Path to reference genome FASTA file (.fasta, .fa, optionally .gz)
    val_transcriptome // string: Path to transcriptome FASTA file for pseudo-alignment (.fasta, .fa, optionally .gz)
    val_gtf // string: Path to gene annotation GTF file (.gtf, optionally .gz)
    val_gtf_isoform // string: Path to extended GTF file for isoform analysis (.gtf, optionally .gz)
    val_bed // string: Path to BED annotation file (.bed)
    val_contaminant_fasta // string: Path to contaminant genome FASTA file for BBSplit (.fasta, .fa, optionally .gz)
    val_rrna_db // string: Path to rRNA database file (.tar.gz) for SortMeRNA
    val_rrna_db_type // string: Runtime mode for SortMeRNA: 'default', 'fast', or 'sensitive'
    val_bbsplit_index // string: Path to BBSplit index folder or tar archive (.tar, .tar.gz)
    val_sortmerna_index // string: Path to SortMeRNA index folder or tar archive (.tar, .tar.gz)
    val_hisat2_index // string: Path to HISAT2 index folder or tar archive (.tar, .tar.gz)
    val_kallisto_index // string: Path to kallisto index file (.idx) - note: kallisto index is a single file, not a folder
    val_salmon_index // string: Path to salmon index folder or tar archive (.tar, .tar.gz)

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

    // Convert string paths to file channels (handle null values)
    ch_fasta = val_fasta ? channel.value(file(val_fasta, checkIfExists: true)) : channel.empty()
    ch_transcriptome = val_transcriptome ? channel.value(file(val_transcriptome, checkIfExists: true)) : channel.empty()
    ch_gtf = val_gtf ? channel.value(file(val_gtf, checkIfExists: true)) : channel.empty()
    ch_gtf_isoform = val_gtf_isoform && "DIU" in params.analysis_method.split(",") ? channel.value(file(val_gtf_isoform, checkIfExists: true)) : channel.empty()
    ch_contaminant_fasta = val_contaminant_fasta && !params.skip_bbsplit ? channel.value(file(val_contaminant_fasta, checkIfExists: true)) : channel.empty()
    ch_rrna_db = val_rrna_db && !params.skip_sortmerna ? channel.value(file(val_rrna_db, checkIfExists: true)) : channel.empty()

    if (params.aligner == "hisat2") {
        // Uncompress genome FASTA if needed
        if (val_fasta.endsWith('.gz')) {
            GUNZIP_GENOME(ch_fasta.map { it -> [[id: it.baseName], it] }.collect())
            ch_fasta_uncompressed = GUNZIP_GENOME.out.gunzip.map { basename, fasta -> fasta }
            ch_versions = ch_versions.mix(GUNZIP_GENOME.out.versions)
        }
        else {
            ch_fasta_uncompressed = ch_fasta
        }
        // Uncompress GTF if needed
        if (val_gtf.endsWith('.gz')) {
            GUNZIP_GTF(ch_gtf.map { it -> [[id: it.baseName], it] }.collect())
            ch_gtf_uncompressed = GUNZIP_GTF.out.gunzip.map { basename, gtf -> gtf }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        }
        else {
            ch_gtf_uncompressed = ch_gtf
        }

        // === HISAT2 INDEX ===
        if (val_hisat2_index) {
            // Use provided index
            if (val_hisat2_index.endsWith('.tar.gz')) {
                UNTAR_HISAT2_INDEX(channel.value(file(val_hisat2_index, checkIfExists: true)))
                ch_hisat2_index = UNTAR_HISAT2_INDEX.out.untar
                ch_versions = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            }
            else {
                ch_hisat2_index = channel.value(file(val_hisat2_index, checkIfExists: true))
            }
        }
        else {
            // Build HISAT2 index from scratch
            // Extract splice sites and exons
            HISAT2_EXTRACTSPLICESITES(ch_gtf_uncompressed)
            ch_splice_sites = HISAT2_EXTRACTSPLICESITES.out.splice_sites
            ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

            HISAT2_EXTRACTEXONS(ch_gtf_uncompressed)
            ch_exon_sites = HISAT2_EXTRACTEXONS.out.exons_sites
            ch_versions = ch_versions.mix(HISAT2_EXTRACTEXONS.out.versions)

            // Build index
            HISAT2_BUILD(ch_fasta_uncompressed, ch_splice_sites, ch_exon_sites)
            ch_hisat2_index = HISAT2_BUILD.out.index
            ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    // Convert GTF to BED for RSeQC
    if (params.rseqc_modules && params.aligner == "hisat2") {
        if (val_bed) {
            ch_bed = channel.value(file(val_bed, checkIfExists: true))
        }
        else {
            GXF2BED(ch_gtf)
            ch_bed = GXF2BED.out.bed
            ch_versions = ch_versions.mix(GXF2BED.out.versions)
        }
    }

    // === SALMON INDEX ===
    if (val_salmon_index) {
        // Use provided index
        if (val_salmon_index.endsWith('.tar.gz')) {
            UNTAR_SALMON_INDEX(channel.value(file(val_salmon_index, checkIfExists: true)))
            ch_salmon_index = UNTAR_SALMON_INDEX.out.untar
            ch_versions = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        }
        else {
            ch_salmon_index = channel.value(file(val_salmon_index, checkIfExists: true))
        }
    }
    else {
        // Build salmon index from transcriptome
        SALMON_INDEX(ch_transcriptome)
        ch_salmon_index = SALMON_INDEX.out.index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    }

    // === KALLISTO INDEX ===
    if (val_kallisto_index) {
        // Use provided index
        if (val_kallisto_index.endsWith('.gz')) {
            GUNZIP_KALLISTO_INDEX(channel.value(file(val_kallisto_index, checkIfExists: true)))
            ch_kallisto_index = GUNZIP_KALLISTO_INDEX.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_KALLISTO_INDEX.out.versions)
        }
        else {
            ch_kallisto_index = channel.value(file(val_kallisto_index, checkIfExists: true))
        }
    }
    else if (params.pseudo_aligner == "kallisto") {
        // Build kallisto index from transcriptome
        KALLISTO_INDEX(ch_transcriptome)
        ch_kallisto_index = KALLISTO_INDEX.out.index
        ch_versions = ch_versions.mix(KALLISTO_INDEX.out.versions)
    }

    // === SORTMERNA INDEX ===
    if (!params.skip_sortmerna) {
        if (val_sortmerna_index) {
            // Use provided index
            if (val_sortmerna_index.endsWith('.tar.gz')) {
                UNTAR_SORTMERNA_INDEX(channel.value(file(val_sortmerna_index, checkIfExists: true)))
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX.out.untar
                ch_versions = ch_versions.mix(UNTAR_SORTMERNA_INDEX.out.versions)
            }
            else {
                ch_sortmerna_index = channel.value(file(val_sortmerna_index, checkIfExists: true))
            }
        }
        else {
            // Build sortmerna index from rRNA database
            UNTAR_SORTMERNA_INDEX(ch_rrna_db)
            ch_rrna_db_fasta = UNTAR_SORTMERNA_INDEX.out.untar
            ch_versions = ch_versions.mix(UNTAR_SORTMERNA_INDEX.out.versions)

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
    }

    // === BBSPLIT INDEX ===
    if (!params.skip_bbsplit && val_contaminant_fasta) {
        if (val_bbsplit_index) {
            // Use provided index
            if (val_bbsplit_index.endsWith('.tar.gz')) {
                UNTAR_BBSPLIT_INDEX(channel.value(file(val_bbsplit_index, checkIfExists: true)))
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX.out.untar
                ch_versions = ch_versions.mix(UNTAR_BBSPLIT_INDEX.out.versions)
            }
            else {
                ch_bbsplit_index = channel.value(file(val_bbsplit_index, checkIfExists: true))
            }
        }
        else {
            // Build bbsplit index from fasta files
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
