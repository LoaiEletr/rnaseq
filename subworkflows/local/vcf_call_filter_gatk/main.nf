//
// Germline variant calling and refinement for RNA-seq data
//
// Performs: GATK4 HaplotypeCaller, BCFtools merge, Tabix indexing, and GATK4 VariantFiltration
//

include { GATK4_HAPLOTYPECALLER } from '../../../modules/local/gatk4/haplotypecaller'
include { BCFTOOLS_MERGE } from '../../../modules/local/bcftools/merge'
include { GATK4_VARIANTFILTRATION } from '../../../modules/local/gatk4/variantfiltration'
include { TABIX_TABIX } from '../../../modules/local/tabix/tabix'

workflow VCF_CALL_FILTER_GATK {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_bai // channel: [ val(meta), [ bai ] ]
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    ch_fai // channel: [ val(meta), [ fai ] ]
    ch_dict // channel: [ val(meta), [ dict ] ]
    ch_intervals // channel: [ val(meta), [ intervals ] ]
    ch_bed // channel: [ val(meta), [ bed ] ]
    ch_dbsnp // channel: [ val(meta), [ dbsnp ] ]
    ch_dbsnp_tbi // channel: [ val(meta), [ dbsnp_tbi ] ]
    skip_variantfiltration // boolean: true/false

    main:

    // 1. Call variants per-interval/per-sample using HaplotypeCaller
    GATK4_HAPLOTYPECALLER(
        ch_bam.join(ch_bai).combine(ch_intervals.map { it[1] }),
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_tbi,
    )

    // Reshape channels to group by condition for merging
    def vcfs_per_condition = GATK4_HAPLOTYPECALLER.out.vcf.map { meta, vcf -> tuple(meta + [id: meta.condition], vcf) }.groupTuple()
    def vcfs_tbi_per_condition = GATK4_HAPLOTYPECALLER.out.tbi.map { meta, tbi -> tuple(meta + [id: meta.condition], tbi) }.groupTuple()

    // 2. Merge scattered VCFs into a single file per condition
    BCFTOOLS_MERGE(
        vcfs_per_condition,
        vcfs_tbi_per_condition,
        ch_bed,
    )

    // 3. Index raw VCF files
    TABIX_TABIX(
        BCFTOOLS_MERGE.out.vcf
    )

    // 4. Apply hard filters to the raw variants
    ch_vcf_final = BCFTOOLS_MERGE.out.vcf
    ch_tbi_final = TABIX_TABIX.out.index

    if (!skip_variantfiltration) {
        GATK4_VARIANTFILTRATION(
            BCFTOOLS_MERGE.out.vcf.join(TABIX_TABIX.out.index),
            ch_fasta,
            ch_fai,
            ch_dict,
        )
        ch_vcf_final = GATK4_VARIANTFILTRATION.out.vcf
    }

    emit:
    vcf = ch_vcf_final // channel: [ val(meta), [ vcf ] ]
    tbi = ch_tbi_final // channel: [ val(meta), [ tbi ] ]
}
