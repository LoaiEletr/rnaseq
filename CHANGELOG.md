# LoaiEletr/rnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - 2026-03-15

Initial release of LoaiEletr/rnaseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

#### Core Pipeline Infrastructure

- Complete end-to-end RNA-seq analysis pipeline orchestration with modular architecture
- Comprehensive genome configuration system with species-specific parameters
- Multi-sample samplesheet processing with 6-column format support
- RO-Crate metadata generation for workflow provenance
- nf-test suite with comprehensive pipeline testing and snapshot validation
- CI/CD configuration with automated testing

#### Preprocessing & Quality Control

- **FastQC**: Initial quality assessment module
- **Cutadapt**: Adapter trimming with quality trimming and sequencer-specific configuration
- **BBMap BBSplit**: Host-contaminant read separation module
- **SortMeRNA**: rRNA removal module
- **UMI-tools suite**: UMI extraction and deduplication modules
- **SeqKit suite**: Sequence sampling and statistics modules
- **Gunzip**: File decompression utility module
- **MultiQC**: Comprehensive pipeline quality reporting with custom configuration support

#### Alignment & Quantification

- **HISAT2 suite**: Complete RNA-Seq alignment pipeline including:
  - Extract exons and splice sites from GTF
  - HISAT2 genome indexing
  - HISAT2 alignment with read group tagging
- **Salmon suite**: Transcript quantification with indexing module
- **Kallisto suite**: Transcript quantification with indexing module
- **Subread featureCounts**: Gene-level counting module
- **Samtools suite**: BAM sorting and indexing modules
- **make_bamlist**: BAM file management module
- **Tximport**: Transcript abundance import module
- **Merge counts**: Combined quantification matrix generation
- **Unified counts matrix generation subworkflow**

#### Quality Assessment (RSeQC Suite)

- **bamstat**: BAM statistics analysis
- **genebodycoverage**: Gene body coverage analysis with heatmap output
- **innerdistance**: Fragment size distribution analysis
- **junctionannotation**: Splice junction analysis with log file output
- **readdistribution**: Genomic feature distribution analysis
- **readduplication**: Duplication rate analysis
- **TIN**: Transcript integrity number calculation
- **inferexperiment**: Strandness inference
- **Comprehensive RSeQC quality control subworkflow**

#### Strandness Detection

- RNA-Seq strandness identification subworkflow
- Library type inference from alignment logs

#### Differential Expression Analysis

- **DESeq2**: Comprehensive differential expression module
- **edgeR-voom**: Differential expression with voom transformation
- **limma_de**: Linear model-based differential expression
- **maSigPro**: Time-course differential expression analysis
- **tx2gene**: Gene identifier mapping module
- **Differential expression visualization module**

#### Functional Enrichment Analysis

- **GO enrichment**: Comprehensive Gene Ontology analysis
- **KEGG pathway**: Pathway enrichment analysis
- **MSigDB GSEA**: Gene set enrichment analysis
- **STRING PPI**: Protein-protein interaction network analysis
- **Comprehensive enrichment analysis suite**

#### Co-expression Network Analysis

- **WGCNA**: Weighted gene co-expression network analysis with sft_r2_threshold parameter
- **Network analysis and enrichment subworkflow**

#### Alternative Splicing Analysis

- **rMATS suite**: Complete alternative splicing analysis including:
  - rMATS for splicing event detection
  - filter_sigevents for result filtering
  - rmats2sashimiplot for splicing visualization
- **IsoformSwitchAnalyzeR**: Isoform switching and alternative splicing analysis
- **DEXSeq suite**: Differential exon usage analysis including:
  - DEXSeq prepare annotation
  - DEXSeq count module
  - DEXSeq DEU module
- **Isoform switch analysis enrichment subworkflow**
- **Differential splicing analysis subworkflow**

#### Variant Calling Pipeline (GVC)

- **GATK4 suite**: Complete variant calling modules including:
  - CreateSequenceDictionary
  - IntervalListTools
  - SplitNCigarReads (RNA-seq specific)
  - BaseRecalibrator
  - ApplyBQSR
  - HaplotypeCaller
  - VariantFiltration
- **bcftools suite**: Variant processing including:
  - merge for VCF concatenation
  - index for VCF indexing
  - isec for variant intersection
- **SnpEff suite**: Variant annotation including:
  - SnpEff download for database installation
  - SnpEff annotation for variant effect prediction
- **Tabix**: Genomic file indexing module
- **Picard MarkDuplicates**: BAM deduplication module
- **BAM deduplication and recalibration subworkflow**
- **Variant calling and filtering subworkflow**
- **Comprehensive Germline Variant Calling pipeline integration**

#### File Format Conversion & Utilities

- **bedtools-sort**: BED file sorting module
- **gtf2bed**: GTF to BED conversion module
- **untar**: Archive extraction utility
- **Conditional execution guards** for error prevention

### `Fixed`

#### Module & Workflow Fixes

- Read group tagging in HISAT2 alignment BAM
- Single-end strandness codes for Kallisto and Salmon quantification
- Unstranded library type detection in HISAT2
- Version capture in versions.yml outputs across all modules
- Stub block outputs for single and paired-end modules
- Conditional execution guards to prevent null errors
- Remove def args from module definitions
- Channel syntax correction (Channel.empty() to channel.empty())
- tx2gene gene identifier correction (ensembl_gene_id vs external_gene_name)
- GTF channel type correction in BAM processing workflow
- Count matrix detection for different input types in edgeR-voom
- Volatile outputs removed from test snapshots (VSD objects, CSV files)
- Strandness output removal from RSeQC inferexperiment
- Proper version collection with .first() in strandness workflow
- rMATS integration issues resolved
- Null check for VST data in visualization script
- MultiQC compatibility with HISAT2 summary format
- Salmon strandness detection log capture for reliable library type inference
- dbsnp parameter made optional in HaplotypeCaller
- Multiple known sites file parsing correction
- Input filename conflicts resolved with stageAs in samtools

### `Changed`

#### Core Infrastructure

- Repository branding and ownership standardization to LoaiEletr
- Pipeline manifest renamed to Loai3tr
- Nextflow version updated to 25.10.2
- nf-test plugin updated from 0.0.3 to 0.0.9
- Template updates from nf-core/tools version 3.5.2
- nf-core-utils plugin added for version tracking

#### Module Standardization

- Module directories renamed to follow nf-core naming conventions
- All modules standardized to nf-core best practices
- Version capturing and I/O formats standardized across all modules
- GATK4 and Picard modules refactored to match best practices
- Module testing standardized with ext.args config files

#### Configuration

- Pipeline configuration and validation streamlined
- Subworkflow outputs consolidated in modules.config
- Default publishDir disabled for better control
- Time allocation increased for low-resource processes
- Tool-specific arguments added to modules.config
- GVC and UMI processing parameters added to config and schema
- GVC-specific reference files added to schema and config

#### Subworkflow Enhancements

- Preprocessing workflow replaced with BBSplit integration
- prepare_genome workflow expanded with comprehensive index building
- GVC parameters integrated into prepare_genome workflow
- Recalibration table output added to BAM dedup subworkflow
- Interval splitting control added to BAM recalibration workflow
- Separate tabix indexing added for filtered VCFs
- HISAT2 subworkflow (deprecated) removed
- Salmon index building enhanced with pre-built indices support
- DEXSeq GFF preparation consolidated in genome setup

#### Test Improvements

- MD5-sensitive files excluded from snapshots
- Inconsistent MD5 snapshots fixed across modules
- DESeq2 test snapshots stabilized by excluding volatile outputs
- Maximum shard count increased for nf-test
- Test data repository path updated

#### Documentation & Metadata

- Tool references and citations updated for QC modules
- GATK, bcftools, and Picard citations added
- Nextflow version badge synced with config
- Parameter descriptions improved with MSigDB category
- RO-Crate metadata updated with correct repository name
- Sample sheet example updated to match 6-column format
- Module documentation enhanced with meta.yml files

### `Dependencies`

#### Core Tools

- Nextflow >=25.10.0
- nf-core/tools 3.5.2 template

#### Analysis Tools Added

- FastQC v0.12.1
- Cutadapt v4.4
- BBMap v39.01
- SortMeRNA v4.3.6
- UMI-tools v1.1.2
- HISAT2 v2.2.1
- Salmon v1.10.0
- Kallisto v0.48.0
- Subread v2.0.6
- Samtools v1.19.2
- RSeQC v5.0.1
- DESeq2 v1.42.0
- edgeR v4.0.0
- limma v3.58.0
- maSigPro v1.74.0
- WGCNA v1.72
- STRINGdb v2.14.0
- rMATS v4.1.2
- rmats2sashimiplot v2.0.5
- IsoformSwitchAnalyzeR v2.2.0
- DEXSeq v1.48.0
- GATK4 v4.5.0.0
- bcftools v1.19
- SnpEff v5.2
- Picard v3.1.1
- Tabix v1.19
- MultiQC v1.21
- SeqKit v2.5.1
- Bedtools v2.31.0
- R v4.3.0 with comprehensive Bioconductor packages

### `Deprecated`

- Old HISAT2 subworkflow replaced with standardized version
