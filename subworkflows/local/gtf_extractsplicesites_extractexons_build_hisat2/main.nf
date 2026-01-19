//
// Extract splice sites and exons from GTF for HISAT2 indexing
//

include { HISAT2_EXTRACTEXONS } from '../../../modules/local/hisat2/extractexons/main.nf'
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/local/hisat2/extractsplicesites/main.nf'
include { HISAT2_BUILD } from '../../../modules/local/hisat2/build/main.nf'

workflow GTF_EXTRACTSPLICESITES_EXTRACTEXONS_BUILD_HISAT2 {
    take:
    ch_gtf // channel: [ gtf ]
    ch_ref // channel: [ ref ]

    main:

    ch_versions = Channel.empty()

    // Extract splice sites from GTF
    HISAT2_EXTRACTSPLICESITES(ch_gtf)
    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

    // Extract exons from GTF
    HISAT2_EXTRACTEXONS(ch_gtf)
    ch_versions = ch_versions.mix(HISAT2_EXTRACTEXONS.out.versions)

    // Build HISAT2 index with reference genome, splice sites, and exons
    HISAT2_BUILD(
        ch_ref,
        HISAT2_EXTRACTSPLICESITES.out.splice_sites,
        HISAT2_EXTRACTEXONS.out.exons_sites,
    )
    ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)

    emit:
    splicesites = HISAT2_EXTRACTSPLICESITES.out.splice_sites // channel: [ splice_sites ]
    exonssites = HISAT2_EXTRACTEXONS.out.exons_sites // channel: [ exons_sites ]
    index = HISAT2_BUILD.out.index // channel: [ index ]
    versions = ch_versions // channel: [ versions.yml ]
}
