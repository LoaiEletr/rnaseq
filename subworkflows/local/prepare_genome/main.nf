//
// Prepare reference genome and build alignment indices for RNA-seq and Germline analysis
// Processes: Genome/GTF decompression, GTF to BED conversion (RSeQC), GTF to GFF conversion (DEXSeq),
// and comprehensive index management for multiple analysis tools.
// Builds indices for: HISAT2 (genome alignment), Kallisto/Salmon (pseudo-alignment),
// SortMeRNA (rRNA filtering), BBMap bbsplit (contaminant removal), and GATK4 (Variant Calling).
// Supports both building indices from scratch or using pre-built indices from compressed archives.
//

include {
    GUNZIP as GUNZIP_GENOME ;
    GUNZIP as GUNZIP_GTF ;
    GUNZIP as GUNZIP_GFF ;
    GUNZIP as GUNZIP_BED ;
    GUNZIP as GUNZIP_KALLISTO_INDEX ;
    GUNZIP as GUNZIP_FAI ;
    GUNZIP as GUNZIP_DICT ;
    GUNZIP as GUNZIP_INTERVALLIST
} from '../../../modules/local/gunzip'
include {
    UNTAR as UNTAR_BBSPLIT_INDEX ;
    UNTAR as UNTAR_SORTMERNA_INDEX ;
    UNTAR as UNTAR_HISAT2_INDEX ;
    UNTAR as UNTAR_SALMON_INDEX ;
    UNTAR as UNTAR_SNPEFF_DB
} from '../../../modules/local/untar'

include {
    TABIX_TABIX as TABIX_TABIX_KNOWNSITES ;
    TABIX_TABIX as TABIX_TABIX_DBSNP
} from '../../../modules/local/tabix/tabix'

include { GXF2BED } from '../../../modules/local/gxf2bed'
include { KALLISTO_INDEX } from '../../../modules/local/kallisto/index'
include { SALMON_INDEX } from '../../../modules/local/salmon/index'
include { HISAT2_EXTRACTEXONS } from '../../../modules/local/hisat2/extractexons'
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/local/hisat2/extractsplicesites'
include { HISAT2_BUILD } from '../../../modules/local/hisat2/build'
include { SORTMERNA as SORTMERNA_INDEX } from '../../../modules/local/sortmerna'
include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_INDEX } from '../../../modules/local/bbmap/bbsplit'
include { DEXSEQ_PREPAREANNOTATION } from '../../../modules/local/dexseq/prepareannotation'
include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools/faidx'
include { SNPEFF_DOWNLOAD } from '../../../modules/local/snpeff/download'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/local/gatk4/createsequencedictionary'
include { GATK4_BEDTOINTERVALLIST } from '../../../modules/local/gatk4/bedtointervallist'
include { GATK4_INTERVALLISTTOOLS } from '../../../modules/local/gatk4/intervallisttools'

workflow PREPARE_GENOME {
    take:
    val_fasta // string: Path to reference genome FASTA file
    val_transcriptome // string: Path to transcriptome FASTA file
    val_gtf // string: Path to gene annotation GTF file
    val_gtf_isoform // string: Path to extended GTF file for isoform analysis
    val_bed // string: Path to BED annotation file
    val_gff // string: Path to GFF annotation file
    val_contaminant_fasta // string: Path to contaminant genome FASTA file
    val_rrna_db // string: Path to rRNA database file (.tar.gz)
    val_rrna_db_type // string: Runtime mode for SortMeRNA: 'default', 'fast', or 'sensitive'
    val_aggregation // boolean: aggregation method for exon counting (DEXSeq)
    val_bbsplit_index // string: Path to BBSplit index folder or tar archive
    val_sortmerna_index // string: Path to SortMeRNA index folder or tar archive
    val_hisat2_index // string: Path to HISAT2 index folder or tar archive
    val_kallisto_index // string: Path to kallisto index file (.idx)
    val_salmon_index // string: Path to salmon index folder or tar archive
    val_fai // string: Path to Samtools Fasta Index (.fai)
    val_dict // string: Path to GATK Sequence Dictionary (.dict)
    val_intervallist // string: Path to GATK Interval List (.interval_list)
    val_dbsnp // string: Path to dbSNP database file
    val_dbsnp_tbi // string: Path to dbSNP database index file
    val_knownsites // string: Path to known sites file
    val_knownsites_tbi // string: Path to known sites index file
    val_snpeff_db // string: Path to SnpEff database folder or tar archive
    val_snpeff_genome // string: SnpEff genome build version (e.g., 'GRCh38.105')

    main:
    ch_fasta_uncompressed = channel.empty()
    ch_gtf_uncompressed = channel.empty()
    ch_rrna_db_fasta = channel.empty()
    ch_hisat2_index = channel.empty()
    ch_kallisto_index = channel.empty()
    ch_salmon_index = channel.empty()
    ch_sortmerna_index = channel.empty()
    ch_bbsplit_index = channel.empty()
    ch_bed = channel.empty()
    ch_gff = channel.empty()
    ch_fai = channel.empty()
    ch_dict = channel.empty()
    ch_intervallist = channel.empty()
    ch_intervals_split = channel.empty()
    ch_snpeff_db = channel.empty()
    ch_dbsnp_tbi = channel.empty()
    ch_knownsites_tbi = channel.empty()

    // Convert string paths to file channels [meta, file]
    ch_fasta = val_fasta ? channel.value(file(val_fasta, checkIfExists: true)).map { [[id: it.baseName], it] }.collect() : channel.empty()
    ch_transcriptome = val_transcriptome ? channel.value(file(val_transcriptome, checkIfExists: true)).map { [[id: it.baseName], it] }.collect() : channel.empty()
    ch_gtf = val_gtf ? channel.value(file(val_gtf, checkIfExists: true)).map { [[id: it.baseName], it] }.collect() : channel.empty()
    ch_gtf_isoform = val_gtf_isoform && ("DIU" in params.analysis_method.split(",") || "AS" in params.analysis_method.split(",")) && params.pseudo_aligner in ["kallisto", "salmon"] ? channel.value(file(val_gtf_isoform, checkIfExists: true)).map { [[id: it.baseName], it] }.collect() : channel.empty()
    ch_contaminant_fasta = val_contaminant_fasta && !params.skip_bbsplit ? channel.value(file(val_contaminant_fasta, checkIfExists: true)).map { [[id: it.baseName], it] }.collect() : channel.empty()
    ch_rrna_db = val_rrna_db && !params.skip_sortmerna ? channel.value(file(val_rrna_db, checkIfExists: true)).map { [[id: it.baseName], it] }.collect() : channel.empty()
    ch_dbsnp = val_dbsnp && "GVC" in params.analysis_method.split(",") ? channel.value(file(val_dbsnp, checkIfExists: true)).map { [[id: it.baseName], it] }.collect() : channel.empty()
    ch_knownsites = val_knownsites && "GVC" in params.analysis_method.split(",") ? channel.fromPath(val_knownsites.split(','), checkIfExists: true).collect().map { files -> [[id: 'knownsites'], files] } : channel.empty()
    ch_snpeff_genome = val_snpeff_genome && "GVC" in params.analysis_method.split(",") ? channel.value([[id: "${val_snpeff_genome}"], val_snpeff_genome]) : channel.empty()

    if (params.aligner == "hisat2") {
        // Uncompress genome FASTA
        if (val_fasta.endsWith('.gz')) {
            GUNZIP_GENOME(ch_fasta)
            ch_fasta_uncompressed = GUNZIP_GENOME.out.gunzip
        }
        else {
            ch_fasta_uncompressed = ch_fasta
        }
        // Uncompress GTF
        if (val_gtf.endsWith('.gz')) {
            GUNZIP_GTF(ch_gtf)
            ch_gtf_uncompressed = GUNZIP_GTF.out.gunzip
        }
        else {
            ch_gtf_uncompressed = ch_gtf
        }

        // === HISAT2 INDEX ===
        if (val_hisat2_index) {
            if (val_hisat2_index.endsWith('.tar.gz')) {
                UNTAR_HISAT2_INDEX(channel.value(file(val_hisat2_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                ch_hisat2_index = UNTAR_HISAT2_INDEX.out.untar
            }
            else {
                ch_hisat2_index = channel.value(file(val_hisat2_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
            }
        }
        else {
            HISAT2_EXTRACTSPLICESITES(ch_gtf_uncompressed)
            HISAT2_EXTRACTEXONS(ch_gtf_uncompressed)
            HISAT2_BUILD(ch_fasta_uncompressed, HISAT2_EXTRACTSPLICESITES.out.splice_sites, HISAT2_EXTRACTEXONS.out.exons_sites)
            ch_hisat2_index = HISAT2_BUILD.out.index
        }

        // === BED ANNOTATION ===
        if ((params.rseqc_modules ? params.rseqc_modules.split(",").any { it in ["bam_stat", "genebody_coverage", "infer_experiment", "inner_distance", "junction_annotation", "read_distribution", "read_duplication", "tin"] } : null) || "GVC" in params.analysis_method.split(",")) {
            if (val_bed) {
                if (val_bed.endsWith('.gz')) {
                    GUNZIP_BED(channel.value(file(val_bed, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                    ch_bed = GUNZIP_BED.out.gunzip
                }
                else {
                    ch_bed = channel.value(file(val_bed, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
                }
            }
            else {
                GXF2BED(ch_gtf)
                ch_bed = GXF2BED.out.bed
            }
        }

        // === GFF ANNOTATION (DEXSeq) ===
        if ("DEU" in params.analysis_method.split(",")) {
            if (val_gff) {
                if (val_gff.endsWith('.gz')) {
                    GUNZIP_GFF(channel.value(file(val_gff, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                    ch_gff = GUNZIP_GFF.out.gunzip
                }
                else {
                    ch_gff = channel.value(file(val_gff, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
                }
            }
            else {
                DEXSEQ_PREPAREANNOTATION(ch_gtf, val_aggregation)
                ch_gff = DEXSEQ_PREPAREANNOTATION.out.gff
            }
        }

        // === GERMLINE VARIANT CALLING (GVC) PREP ===
        if ("GVC" in params.analysis_method.split(",")) {
            // FASTA INDEX
            if (val_fai) {
                if (val_fai.endsWith('.gz')) {
                    GUNZIP_FAI(channel.value(file(val_fai, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                    ch_fai = GUNZIP_FAI.out.gunzip
                }
                else {
                    ch_fai = channel.value(file(val_fai, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
                }
            }
            else {
                SAMTOOLS_FAIDX(ch_fasta_uncompressed)
                ch_fai = SAMTOOLS_FAIDX.out.fai
            }

            // DICTIONARY
            if (val_dict) {
                if (val_dict.endsWith('.gz')) {
                    GUNZIP_DICT(channel.value(file(val_dict, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                    ch_dict = GUNZIP_DICT.out.gunzip
                }
                else {
                    ch_dict = channel.value(file(val_dict, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
                }
            }
            else {
                GATK4_CREATESEQUENCEDICTIONARY(ch_fasta_uncompressed)
                ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
            }

            // INTERVAL LIST
            if (val_intervallist) {
                if (val_intervallist.endsWith('.gz')) {
                    GUNZIP_INTERVALLIST(channel.value(file(val_intervallist, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                    ch_intervallist = GUNZIP_INTERVALLIST.out.gunzip
                }
                else {
                    ch_intervallist = channel.value(file(val_intervallist, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
                }
            }
            else {
                GATK4_BEDTOINTERVALLIST(ch_bed, ch_dict)
                ch_intervallist = GATK4_BEDTOINTERVALLIST.out.interval_list
            }

            // SPLITTING
            ch_intervals_split = ch_intervallist
            if (!params.skip_interval_splitting) {
                GATK4_INTERVALLISTTOOLS(ch_intervallist)
                ch_intervals_split = GATK4_INTERVALLISTTOOLS.out.interval_list
            }

            // SNPEFF DB
            if (val_snpeff_db) {
                if (val_snpeff_db.endsWith('.tar.gz')) {
                    UNTAR_SNPEFF_DB(channel.value(file(val_snpeff_db, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                    ch_snpeff_db = UNTAR_SNPEFF_DB.out.untar
                }
                else {
                    ch_snpeff_db = channel.value(file(val_snpeff_db, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
                }
            }
            else {
                SNPEFF_DOWNLOAD(
                    ch_snpeff_genome
                )
                ch_snpeff_db = SNPEFF_DOWNLOAD.out.cache
            }

            // DBSNP INDEX
            if (val_dbsnp) {

                if (val_dbsnp_tbi) {
                    ch_dbsnp_tbi = channel.value(file(val_dbsnp_tbi, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
                }
                else {
                    TABIX_TABIX_DBSNP(
                        ch_dbsnp
                    )
                    ch_dbsnp_tbi = TABIX_TABIX_DBSNP.out.index
                }
            }

            // KNOWNSITES INDEX
            if (val_knownsites) {

                if (val_knownsites_tbi) {
                    ch_knownsites_tbi = channel.fromPath(val_knownsites_tbi.split(','), checkIfExists: true).collect().map { indices -> [[id: 'knownsites'], indices] }
                }
                else {
                    TABIX_TABIX_KNOWNSITES(
                        ch_knownsites
                    )
                    ch_knownsites_tbi = TABIX_TABIX_KNOWNSITES.out.index.collect { it[1] }.map { indices -> [[id: 'knownsites'], indices] }
                }
            }
        }
    }

    // === SALMON INDEX ===
    if (val_salmon_index) {
        if (val_salmon_index.endsWith('.tar.gz')) {
            UNTAR_SALMON_INDEX(channel.value(file(val_salmon_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
            ch_salmon_index = UNTAR_SALMON_INDEX.out.untar
        }
        else {
            ch_salmon_index = channel.value(file(val_salmon_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
        }
    }
    else {
        SALMON_INDEX(ch_transcriptome)
        ch_salmon_index = SALMON_INDEX.out.index
    }

    // === KALLISTO INDEX ===
    if (val_kallisto_index) {
        if (val_kallisto_index.endsWith('.gz')) {
            GUNZIP_KALLISTO_INDEX(channel.value(file(val_kallisto_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
            ch_kallisto_index = GUNZIP_KALLISTO_INDEX.out.gunzip
        }
        else {
            ch_kallisto_index = channel.value(file(val_kallisto_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
        }
    }
    else if (params.pseudo_aligner == "kallisto") {
        KALLISTO_INDEX(ch_transcriptome)
        ch_kallisto_index = KALLISTO_INDEX.out.index
    }

    // === SORTMERNA INDEX ===
    if (!params.skip_sortmerna) {
        if (val_sortmerna_index) {
            if (val_sortmerna_index.endsWith('.tar.gz')) {
                UNTAR_SORTMERNA_INDEX(channel.value(file(val_sortmerna_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX.out.untar
            }
            else {
                ch_sortmerna_index = channel.value(file(val_sortmerna_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
            }
        }
        else {
            UNTAR_SORTMERNA_INDEX(ch_rrna_db)
            ch_rrna_db_fasta = UNTAR_SORTMERNA_INDEX.out.untar
            SORTMERNA_INDEX([[id: 'index'], []], [[:], []], ch_rrna_db_fasta, val_rrna_db_type, true)
            ch_sortmerna_index = SORTMERNA_INDEX.out.index
        }
    }

    // === BBSPLIT INDEX ===
    if (!params.skip_bbsplit && val_contaminant_fasta) {
        if (val_bbsplit_index) {
            if (val_bbsplit_index.endsWith('.tar.gz')) {
                UNTAR_BBSPLIT_INDEX(channel.value(file(val_bbsplit_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect())
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX.out.untar
            }
            else {
                ch_bbsplit_index = channel.value(file(val_bbsplit_index, checkIfExists: true)).map { [[id: it.baseName], it] }.collect()
            }
        }
        else {
            BBMAP_BBSPLIT_INDEX([[id: 'index'], []], [[:], []], ch_fasta, ch_contaminant_fasta, true)
            ch_bbsplit_index = BBMAP_BBSPLIT_INDEX.out.index
        }
    }

    emit:
    fasta_compressed = ch_fasta // channel: [ val(meta), [ fasta ] ]
    fasta_uncompressed = ch_fasta_uncompressed // channel: [ val(meta), [ fasta ] ]
    transcriptome = ch_transcriptome // channel: [ val(meta), [ fasta ] ]
    gtf_compressed = ch_gtf // channel: [ val(meta), [ gtf ] ]
    gtf_uncompressed = ch_gtf_uncompressed // channel: [ val(meta), [ gtf ] ]
    gtf_isoform = ch_gtf_isoform // channel: [ val(meta), [ gtf ] ]
    contaminant_fasta = ch_contaminant_fasta // channel: [ val(meta), [ fasta ] ]
    rrna_db_fasta = ch_rrna_db_fasta // channel: [ val(meta), [ fasta ] ]
    bed = ch_bed // channel: [ val(meta), [ bed ] ]
    gff = ch_gff // channel: [ val(meta), [ gff ] ]
    fai = ch_fai // channel: [ val(meta), [ fai ] ]
    dict = ch_dict // channel: [ val(meta), [ dict ] ]
    intervallist = ch_intervallist // channel: [ val(meta), [ interval_list ] ]
    intervals_split = ch_intervals_split // channel: [ val(meta), [ interval_list ] ]
    knownsites = ch_knownsites // channel: [ val(meta), [ knownsites ] ]
    knownsites_tbi = ch_knownsites_tbi // channel: [ val(meta), [ knownsites_tbi ] ]
    dbsnp = ch_dbsnp // channel: [ val(meta), [ dbsnp ] ]
    dbsnp_tbi = ch_dbsnp_tbi // channel: [ val(meta), [ dbsnp_tbi ] ]
    snpeff_db = ch_snpeff_db // channel: [ val(meta), [ snpeff_db ] ]
    snpeff_genome = ch_snpeff_genome // channel: [ val(meta), [ snpeff_genome ] ]
    hisat2_index = ch_hisat2_index // channel: [ val(meta), [ hisat2_index ] ]
    kallisto_index = ch_kallisto_index // channel: [ val(meta), [ kallisto_index ] ]
    salmon_index = ch_salmon_index // channel: [ val(meta), [ salmon_index ] ]
    sortmerna_index = ch_sortmerna_index // channel: [ val(meta), [ sortmerna_index ] ]
    bbsplit_index = ch_bbsplit_index // channel: [ val(meta), [ bbsplit_index ] ]
}
