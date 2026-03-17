# LoaiEletr/rnaseq: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**LoaiEletr/rnaseq** is a comprehensive bioinformatics pipeline designed for the end-to-end analysis of transcriptomic data. Built using Nextflow DSL2 and the nf-core framework, it ensures maximum reproducibility and scalability across various compute environments.

The workflow ingests a samplesheet containing FASTQ files and performs extensive pre-processing, including quality control, adapter trimming, contamination removal, and ribosomal RNA depletion.

### Branched Analysis Architecture

A key feature of this pipeline is its **branched architecture**. Your choice of alignment tool determines the downstream analysis modules available:

- **Genomic Route (HISAT2)**: Best for users requiring genomic coordinates (BAM files). This route supports Differential Gene Expression (DEG), WGCNA, Alternative Splicing (rMATS), Differential Exon Usage (DEXSeq), and Germline Variant Calling (GVC).
- **Pseudoalignment Route (Salmon/Kallisto)**: Optimized for high-speed transcript quantification. This route supports DEG, WGCNA, and Isoform Switching analysis (via IsoformSwitchAnalyzeR).

### Key Features

- **Automated QC**: Integration of `FastQC`, `RSeQC`, and `MultiQC` for comprehensive reporting at every stage.
- **Flexibility**: Support for Unique Molecular Identifiers (UMIs) to mitigate PCR bias and automatic library strandedness detection.
- **Advanced Analysis**: Beyond standard differential expression, the pipeline offers integrated modules for time-course data (`maSigPro`), co-expression networks (`WGCNA`), and protein-protein interaction (`STRING`).

> [!NOTE]
> For a visual overview of the major workflow steps, please refer to the "tube map" diagram in the [main README](README.md).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Technical Replicates

Note that the pipeline **does not** support automatic concatenation of multi-lane or technical replicates. Each row in the samplesheet must represent a unique biological replicate with a unique `sample_id`.

If you have re-sequenced the same sample across multiple lanes or runs to increase depth, you must manually concatenate the FASTQ files (e.g., using the `cat` command) before creating the samplesheet.

**Example of unsupported configuration:**
Do not use the same `sample_id` for multiple rows:

| sample_id | fastq_1          | fastq_2          |
| --------- | ---------------- | ---------------- |
| SAMPLE_A  | run1_R1.fastq.gz | run1_R2.fastq.gz |
| SAMPLE_A  | run2_R1.fastq.gz | run2_R2.fastq.gz |

Instead, concatenate these files beforehand so that each `sample_id` appears only once in your input file.

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 6 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below.

```csv title="samplesheet.csv"
sample_id,fastq_1,fastq_2,condition,lib_type,sequencer
CONTROL_REP1,/path/to/fastq/files/control_rep1_R1.fastq.gz,/path/to/fastq/files/control_rep1_R2.fastq.gz,control,reverse,NextSeq
CONTROL_REP2,/path/to/fastq/files/control_rep2_R1.fastq.gz,/path/to/fastq/files/control_rep2_R2.fastq.gz,control,reverse,NovaSeq
CONTROL_REP3,/path/to/fastq/files/control_rep3_R1.fastq.gz,/path/to/fastq/files/control_rep3_R2.fastq.gz,control,reverse,NovaSeq
TREATED_REP1,/path/to/fastq/files/treated_rep1_R1.fastq.gz,/path/to/fastq/files/treated_rep1_R2.fastq.gz,treated,forward,MiSeq
TREATED_REP2,/path/to/fastq/files/treated_rep2_R1.fastq.gz,/path/to/fastq/files/treated_rep2_R2.fastq.gz,treated,forward,HiSeq
TREATED_REP3,/path/to/fastq/files/treated_rep3_R1.fastq.gz,/path/to/fastq/files/treated_rep3_R2.fastq.gz,treated,forward,HiSeq
```

Each row represents a biological replicate. The pipeline auto-detects whether a sample is single-end (only `fastq_1` provided) or paired-end (both `fastq_1` and `fastq_2` provided).

| Column      | Description                                                                                                                                                                                                                                                                                                                                                    |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample_id` | Custom sample name (no spaces allowed). Unique identifier for each biological replicate.                                                                                                                                                                                                                                                                       |
| `fastq_1`   | Full path to FastQ file for read 1. File must be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                                                                                                                                                                                       |
| `fastq_2`   | Full path to FastQ file for read 2 (paired-end only). Leave empty for single-end data.                                                                                                                                                                                                                                                                         |
| `condition` | Experimental condition. Must contain **exactly two unique values** across all samples, with **at least 3 replicates per condition**. The control condition **must** contain the word "control" (case-insensitive, e.g., "control", "CONTROL", "Control_sample"). The other condition can be any user-defined name (e.g., "treated", "knockout", "stimulated"). |
| `lib_type`  | Library strandedness. Must be one of: `forward`, `reverse`, or `auto`. If set to `auto`, the pipeline will subsample reads and use Salmon to automatically infer strandedness for correct alignment.                                                                                                                                                           |
| `sequencer` | Sequencing platform. Must be one of: `HiSeq`, `MiSeq`, `NovaSeq`, or `NextSeq`. Used for platform-specific trimming parameters (e.g., `--nextseq-trim=20` for NovaSeq/NextSeq data).                                                                                                                                                                           |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

**Important requirements:**

- **Two conditions only**: The pipeline supports exactly two experimental conditions (e.g., control vs treated). Do not include more than two unique values in the `condition` column.
- **Minimum 3 replicates per condition**: Each condition must have at least three biological replicates to ensure statistical power for differential analysis.
- **Control condition naming**: Must contain the word "control" (case-insensitive) to properly identify the reference group for comparisons.
- **No multi-lane/technical replicates**: Each row represents a unique biological replicate. Technical replicates (multiple sequencing runs of the same sample) are not supported.
- **Single or paired-end**: Automatically detected based on whether `fastq_2` is provided or empty.
- **Strandedness auto-detection**: Use `auto` in the `lib_type` column if you are unsure about the library strandedness. The pipeline will infer it automatically via Salmon subsampling.

---

## Analysis Methods

The pipeline operates using a branched architecture where your choice of aligner determines the available downstream analysis modules.

> [!IMPORTANT]
> Users must select **exactly one aligner per run** (`--aligner` or `--pseudo_aligner`). The `--analysis_method` parameter accepts a comma-separated list, allowing you to run **multiple compatible analyses simultaneously**.

### Analysis Route Compatibility

| Feature                | HISAT2 Route                       | Pseudo-aligner Route (Salmon/Kallisto)  |
| :--------------------- | :--------------------------------- | :-------------------------------------- |
| **Command**            | `--aligner hisat2`                 | `--pseudo_aligner salmon` or `kallisto` |
| **Primary Output**     | Genomic Alignments (BAM)           | Transcript Quantifications (SF/TSV)     |
| **Compatible Methods** | `DEG`, `WGCNA`, `AS`, `DEU`, `GVC` | `DEG`, `WGCNA`, `DIU`, `AS`             |

### Method Descriptions

- **DEG**: Differential Gene Expression using `DESeq2`, `limma`, or `maSigPro`.
- **WGCNA**: Weighted Gene Co-expression Network Analysis to identify highly correlated gene modules.
- **AS**: Alternative Splicing analysis (via `rMATS` for HISAT2; via `IsoformSwitchAnalyzeR` for Pseudo-aligners).
- **DEU**: Differential Exon Usage using `DEXSeq` (HISAT2 only).
- **DIU**: Differential Isoform Usage using `IsoformSwitchAnalyzeR` (Pseudo-aligner only).
- **GVC**: Germline Variant Calling using the `GATK4` best practices (HISAT2 only).

### Example Usage

#### Genomic Alignment (HISAT2)

To run differential expression, co-expression networks, alternative splicing, and differential exon usage in a single command:

```bash
--aligner hisat2 --analysis_method DEG,WGCNA,AS,DEU
```

#### Transcriptome Pseudo-alignment (Salmon or Kallisto)

To run differential expression, co-expression networks, and both isoform switching and alternative splicing analysis:

```bash
--pseudo_aligner salmon --analysis_method DEG,WGCNA,DIU,AS
```

OR

```bash
--pseudo_aligner kallisto --analysis_method DEG,WGCNA,DIU,AS
```

---

## UMI Extraction

Unique Molecular Identifiers (UMIs) are short sequences added to each DNA molecule before PCR amplification. They are essential for identifying and removing PCR duplicates, ensuring that each read is counted only once during quantification.

To enable UMI processing in the pipeline, use the `--with_umi` flag.

### UMI Parameters

- `--with_umi`: Set to `true` to enable the UMI extraction and deduplication workflow.
- `--skip_umi_extract`: If your reads already have the UMI sequence appended to the read header (e.g., by a sequencing facility), you can skip the extraction step while still performing deduplication.
- `--lib_kit`: Name of the library preparation kit (e.g., `quantseq`, `corall`, `takara`) to use a preset UMI pattern.
- `--umitools_bc_pattern`: Custom UMI barcode pattern. This takes priority over `--lib_kit` if both are provided.
- `--umitools_bc_pattern2`: Used for the second read in paired-end data if the UMI is split or located on both reads.
- `--umi_separator`: The character used to separate the UMI from the read ID in the header.
- `--umi_extract_method`: The method used to extract the barcode. Usually `string` (fixed position) or `regex` (pattern matching).
- `--umi_post_extract_method`: Specific downstream filtering or processing logic applied after the UMI has been moved to the read header.

### Defining the UMI Pattern

The pipeline needs to know where the UMI is located within your reads. You can define this in two ways:

#### 1. Using a Library Kit Preset

If you are using a standard commercial kit, you can use the `--lib_kit` parameter. The pipeline will automatically apply the correct barcode pattern for you:

| Kit Option | Applied Pattern                         |
| :--------- | :-------------------------------------- |
| `quantseq` | `^(?P<umi_1>.{6})(?P<discard_1>.{4}).*` |
| `corall`   | `^(?P<umi_1>.{12}).*`                   |
| `takara`   | `^(?P<umi_1>.{8})(?P<discard_1>.{6}).*` |

#### 2. Using a Custom Pattern (Priority)

If your kit is not listed above or you have a custom design, use `--umitools_bc_pattern`.

> **Note:** If you provide a custom pattern via `--umitools_bc_pattern`, it will take priority and override any selection made via `--lib_kit`.

### Example Usage

```bash
--with_umi --umitools_bc_pattern 'NNNNNN' --umitools_bc_pattern2 'NNNNNN' --umi_separator '_' --umi_extract_method 'string'
```

---

## Adapter Trimming

The pipeline uses [`Cutadapt`](https://cutadapt.readthedocs.io/) to remove sequencing adapters and low-quality bases from the ends of reads. This step is critical for ensuring accurate alignment and reducing false-positive variant calls.

### Trimming Parameters

- `--adapters`: (Optional) Path to a custom FASTA file containing adapter sequences.
- `--skip_trimming`: Set to `true` to skip the adapter trimming step.

> **Note**
> Enabling `--skip_trimming` also automatically skips the **post-trimming QC (FastQC)**. This reduces redundant processing if you are providing pre-processed reads.

### Automatic Adapter Selection

If you do not provide a custom file via `--adapters`, the pipeline intelligently selects the appropriate universal adapter set based on your input data:

1. **Detection**: The pipeline checks if the samplesheet contains single-end or paired-end reads.
2. **Selection**:
   - For **Single-End (SE)**: Uses the presets in `assets/adapters-SE.fa`.
   - For **Paired-End (PE)**: Uses the presets in `assets/adapters-PE.fa`.
3. **Execution**: Cutadapt is then run with these sequences to clean your raw FASTQ data.

### Example Usage

**Using custom adapters:**

```bash
--adapters '/path/to/my_adapters.fa'
```

---

## Reference Genome Configuration

The pipeline requires a reference genome to align your reads. You can either provide local paths to your files or use one of the built-in species presets.

### Reference Parameters

- `--species`: The target species name (e.g., `human`, `mouse`, `rat`).
- `--genome`: The specific genome build version (e.g., `GRCh38`, `GRCm39`, `GRCr8`).
- `--fasta`: (Optional) Manual path to the genome FASTA file.
- `--gtf`: (Optional) Manual path to the gene annotation GTF file.
- `--bed`: (Optional) Manual path to gene annotation BED file.
- `--transcriptome`: (Optional) Manual path to transcriptome FASTA file.
- `--gtf_isoform`: (Optional) Manual path to extended GTF for isoform analysis.

> [!TIP]
> When using `--species` and `--genome` presets, the pipeline automatically retrieves the correct paths for the genome FASTA, transcriptome, GTF (including `IsoformSwitchAnalyzeR` versions), and dbSNP files.

> [!IMPORTANT]
> **Why specify `--species`?** > Even if you provide manual paths to your FASTA and GTF files, you should still specify the `--species` from the supported list below. This ensures the pipeline uses the correct organism-specific database for Differential Expression (DEG) and Functional Enrichment.

### Supported Species & Genome Versions

If you choose to use the built-in presets, you must provide both the `--species` and the corresponding `--genome` version as shown in the table below:

| Species       | Supported `--genome` Version | Description              |
| :------------ | :--------------------------- | :----------------------- |
| `human`       | `GRCh38`                     | Homo sapiens             |
| `mouse`       | `GRCm39`                     | Mus musculus             |
| `rat`         | `GRCr8`                      | Rattus norvegicus        |
| `monkey`      | `Mmul_10`                    | Macaca mulatta           |
| `chicken`     | `GRCg7b`                     | Gallus gallus            |
| `cow`         | `ARS-UCD2.0`                 | Bos taurus               |
| `pig`         | `Sscrofa11.1`                | Sus scrofa               |
| `dog`         | `ROSCfam1.0`                 | Canis lupus familiaris   |
| `yeast`       | `R64-1-1`                    | Saccharomyces cerevisiae |
| `zebrafish`   | `GRCz11`                     | Danio rerio              |
| `fruitfly`    | `BDGP6`                      | Drosophila melanogaster  |
| `worm`        | `WBcel235`                   | Caenorhabditis elegans   |
| `arabidopsis` | `TAIR10`                     | Arabidopsis thaliana     |

### Example Parameters

**Using a preset for Human:**

```bash
--species human --genome GRCh38
```

**Using a preset for Mouse:**

```bash
--species human --genome GRCh38
```

**Using custom local files (Overrides presets):**

```bash
--fasta '/path/to/genome.fa' --gtf '/path/to/annotation.gtf'
```

---

## Contamination Removal (BBSplit)

In many experimental setups (e.g., human tumor xenografts in mice or cell cultures with potential fungal contamination), it is necessary to filter out reads belonging to a non-target species. The pipeline uses [`BBSplit`](https://jgi.doe.gov/data-and-tools/bbtools/) to achieve this.

### Contamination Parameters

- `--contaminant_species`: Automatically filters out reads mapping to a common lab contaminant. Supported presets:
  - `human`, `mouse`, `rat`, `monkey`, `chicken`, `yeast`, `zebrafish`, `fruitfly`, `worm`.
- `--contaminant_fasta`: Manual path to a FASTA file of a specific contaminant not listed in the presets.
- `--bbsplit_index`: (Optional) Path to a pre-computed BBSplit index directory or `.tar.gz` archive (Must be named `bbsplit`).
- `--skip_bbsplit`: Set to `true` to skip the contamination filtering step entirely, even if a contaminant species or FASTA is provided.

### How it Works

1. **Target vs. Contaminant**: The pipeline simultaneously aligns reads to your **Target** (defined by `--species`) and your **Contaminant** (defined by `--contaminant_species`).
2. **Filtering**: Reads that align better to the contaminant genome are discarded.
3. **Downstream**: Only the "clean" reads that align best to your target species are passed to the rest of the pipeline.

> [!IMPORTANT]
> If providing a manual index, the container (directory or `.tar.gz`) must be named **`bbsplit`**. Inside the index, the primary reference files must be prefixed with `primary` to ensure the pipeline correctly identifies the target reads to keep.

### Example Parameters

**Using a pre-built BBSplit index archive:**

```bash
--bbsplit_index '/path/to/bbsplit.tar.gz'
```

**Study: Human cells in a Mouse xenograft**

```bash
--species human --genome GRCh38 --contaminant_species mouse
```

**Study: Mouse cells with suspected Yeast contamination**

```bash
--species mouse --genome GRCm39 --contaminant_species yeast
```

**Using a custom contaminant (e.g., Mycoplasma)**

```bash
--species human --genome GRCh38 --contaminant_fasta '/path/to/mycoplasma.fa'
```

**Disabling BBSplit (Ignoring contaminant settings):**

```bash
--skip_bbsplit
```

---

## rRNA Removal (SortMeRNA)

The pipeline uses [`SortMeRNA`](https://github.com/biocore/sortmerna) to filter out ribosomal RNA (rRNA) sequences. This step is essential for improving the signal-to-noise ratio in your transcriptomic data by removing highly abundant but non-informative rRNA reads.

### rRNA Parameters

- `--rrna_db_type`: Selects the specific SortMeRNA database profile. Available options:
  - `fast` (**Default**): Balanced for speed and low RAM usage (99.88% accuracy).
  - `default`: Standard sensitivity (99.89% accuracy).
  - `sensitive`: Higher clustering identity; runs ~2x slower.
  - `sensitive_rfam`: Maximum sensitivity including full RFAM seed sequences.
- `--sortmerna_index`: (Optional) Path to a pre-computed SortMeRNA index directory or `.tar.gz` archive (Must be named `index`).
- `--skip_sortmerna`: Set to `true` to skip rRNA removal entirely.

### Database Composition

The pipeline automatically retrieves and manages the SortMeRNA v4.3.4 database, which includes sequences from:

- **SILVA 138/132**: 16S, 18S, 23S, and 28S.
- **RFAM 14.1**: 5S and 5.8S.

| Database Type    | Reference File Used                      | Accuracy | Speed        |
| :--------------- | :--------------------------------------- | :------- | :----------- |
| `fast`           | `smr_v4.3_fast_db.fasta`                 | 99.888%  | ⚡ Very Fast |
| `default`        | `smr_v4.3_default_db.fasta`              | 99.899%  | 🟢 Normal    |
| `sensitive`      | `smr_v4.3_sensitive_db.fasta`            | 99.907%  | 🟠 Slow      |
| `sensitive_rfam` | `smr_v4.3_sensitive_db_rfam_seeds.fasta` | >99.907% | 🔴 Very Slow |

> [!TIP]
> If you are working on a low-resource machine (e.g., 2 cores, 8GB RAM), using a pre-computed index via `--sortmerna_index` is highly recommended to bypass the memory-intensive indexing phase.

### Example Parameters

**Using a pre-built SortMeRNA index archive:**

```bash
--sortmerna_index '/path/to/index.tar.gz'
```

**Using the high-sensitivity database:**

```bash
--rrna_db_type 'sensitive'
```

**Skipping rRNA removal:**

```bash
--skip_sortmerna
```

---

## Strandedness Inference

The pipeline features an automated detection system to identify the library strandedness of your samples. This is essential for accurate gene quantification, especially in regions with overlapping genes on opposite strands.

### How it Works

The inference step is triggered for any sample where the `lib_type` column in your samplesheet is set to `auto`.

1. **Subsampling**: The pipeline extracts a random subset of **1 million reads** (R1/R2 for paired-end, or R1 for single-end).
2. **Pseudo-mapping**: These reads are mapped to the transcriptome using [`Salmon`](https://salmon.readthedocs.io/).
3. **Classification**: Salmon's `libType` auto-detection logic determines the most likely orientation (e.g., `ISR`, `ISF`, `IU`).
4. **Application**: The detected strandedness is automatically applied to all downstream quantification and analysis steps.

### Samplesheet Configuration

To use this feature, simply set the `lib_type` column to `auto`:

| sample_id | fastq_1          | fastq_2          | condition | lib_type | sequencer |
| :-------- | :--------------- | :--------------- | :-------- | :------- | :-------- |
| CONTROL_1 | path/to/R1.fq.gz | path/to/R2.fq.gz | control   | `auto`   | NextSeq   |

> [!TIP]
> If you already know your library type (e.g., `forward` or `reverse`), it is recommended to specify it directly to save the extra processing time required for subsampling and inference.

### Example Parameters

**Run the pipeline allowing auto-detection for specific samples:**

```bash
--samplesheet 'samplesheet.csv' --species human --genome GRCh38
```

---

## Alignment & Quantification

The pipeline supports three distinct routes for processing reads. You can choose between traditional genomic alignment or high-speed transcriptomic pseudoalignment.

### Alignment Parameters

- `--pseudo_aligner`: Selects the pseudoalignment tool. Options: `kallisto` (**Default**) or `salmon`.
- `--aligner`: Selects a traditional splice-aware aligner. Currently supports: `hisat2`.

> [!CAUTION]
> The pipeline cannot run a pseudo-aligner and a genomic aligner simultaneously. If you wish to use `--aligner hisat2`, you **must** set `--pseudo_aligner null` in your parameters to avoid a system conflict.

### Reference Requirements

The pipeline automatically manages these files based on your `--species` and `--genome` selection, but they must be available in the genome definition:

| Tool         | Required Reference Files                    |
| :----------- | :------------------------------------------ |
| **HISAT2**   | `--fasta` (Genome) and `--gtf` (Annotation) |
| **Salmon**   | `--transcriptome` and `--gtf`               |
| **Kallisto** | `--transcriptome` and `--gtf`               |

### Tool-Specific Options

#### Kallisto (Pseudoalignment)

- **Single-End Mode**: Unlike paired-end reads, single-end reads do not provide internal fragment length information. You must provide:
  - `--fragment_length`: Estimated average fragment length (Default: `200`).
  - `--fragment_length_sd`: Estimated standard deviation of fragment length (Default: `30`).
- **Bootstrapping (Optional)**:
  - `--bootstrap_count`: Number of bootstrap samples to perform (Default: `0`). This is used to estimate technical variance in abundance estimates, which is useful for downstream tools like Sleuth.

### Index Management

For each tool, the pipeline can either **build the index automatically** (using your FASTA/GTF) or use a **pre-computed index** provided by you to save time.

| Tool         | Parameter          | Required Format                                                |
| :----------- | :----------------- | :------------------------------------------------------------- |
| **Kallisto** | `--kallisto_index` | A single file: `.idx` or `.idx.gz`.                            |
| **Salmon**   | `--salmon_index`   | A directory containing the index files OR a `.tar.gz` archive. |
| **HISAT2**   | `--hisat2_index`   | A directory containing `.ht2` files OR a `.tar.gz` archive.    |

> [!IMPORTANT]
> If providing a manual HISAT2 index, the internal filenames must use the **base name** of your reference FASTA.
>
> - _Example:_ If your reference is `genome.fasta`, the index files must be `genome.1.ht2`, `genome.2.ht2`, etc.

### Example Parameters

**Using a pre-built Kallisto index archive:**

```bash
--kallisto_index '/path/to/kallisto_index.idx.gz'
```

**Using a pre-built Salmon index archive:**

```bash
--salmon_index '/path/to/salmon_index.tar.gz'
```

**Using a pre-built HISAT2 index archive:**

```bash
--hisat2_index '/path/to/hisat2_index.tar.gz'
```

---

## Post-Alignment Processing (SAMtools)

If you choose genomic alignment via `HISAT2`, the pipeline performs essential post-processing of the raw alignment data to prepare it for visualization and downstream variant calling.

### Processing Steps

Once the alignment is complete, the pipeline automatically uses [`SAMtools`](https://www.htslib.org/) to:

1.  **Coordinate Sorting**: Reorders the mapped reads by their genomic coordinates (Chromosome and Position).
2.  **Indexing**: Generates a `.bam.bai` index file, allowing tools like IGV or GATK to access specific genomic regions rapidly.

> [!NOTE]
> This step is automatically skipped if you are using a pseudo-aligner (`kallisto` or `salmon`), as those tools quantify transcripts directly without producing coordinate-sorted genomic BAM files.

---

## UMI-based Deduplication (UMI-tools)

For libraries utilizing Unique Molecular Identifiers (UMIs) to mitigate PCR duplication bias, the pipeline integrates [`UMI-tools`](https://github.com/CGATOxford/UMI-tools). This workflow ensures that multiple PCR copies of the same original RNA molecule are collapsed into a single representative read.

### UMI Parameters

- `--with_umi`: Set to `true` to enable the UMI extraction and deduplication workflow.
- `--umi_separator`: The character used to separate the UMI from the read ID in the FASTQ header (e.g., `_` or `:`).
- `--umi_post_extract_method`: Defines specific downstream filtering or processing logic applied after the UMI has been moved to the read header (e.g., `unique`, `percentile`, or `adjacency`).

> [!CAUTION]
> **HISAT2 Requirement**: UMI deduplication is currently only supported when using the `--aligner hisat2` route. This process operates on coordinate-sorted BAM files, which are not produced by the Salmon or Kallisto pseudo-alignment routes.

### Workflow Integration

1. **Extraction**: UMIs are removed from the sequenced read and appended to the read name in the FASTQ file.
2. **Alignment**: Reads are aligned via **HISAT2** using the UMI-tagged headers.
3. **Deduplication**: After alignment, `UMI-tools dedup` uses the mapping coordinates and the UMI sequence to identify and remove PCR duplicates.

### Example Parameters

**Standard UMI processing:**

```bash
--with_umi --umi_separator '_' --umi_post_extract_method 'directional'
```

---

## Matrix Assembly

Matrix assembly is the process of combining individual sample results into a unified data frame. This step is automatically triggered if the `--analysis_method` parameter is provided.

The pipeline selects the assembly route based on your chosen aligner:

### 1. Genomic Route (HISAT2 + featureCounts)

When using `--aligner hisat2`, the pipeline uses `featureCounts` to assign reads to genomic features.

- **Mechanism**: A custom R script aggregates individual `featureCounts` outputs.
- **Output**: A single **Gene Count Matrix** containing raw integer counts.

### 2. Pseudoalignment Route (Salmon / Kallisto)

When using `--pseudo_aligner`, the pipeline utilizes the [`tximport`](https://bioconductor.org/packages/release/bioc/html/tximport.html) R package.

- **Mechanism**: `tximport` summarizes transcript-level abundances (TPM) and estimated counts.
- **Benefit**: Automatically accounts for gene-length changes across samples, improving downstream accuracy.

### Assembly & Analysis Parameters

- `--analysis_method`: Defines the downstream statistical processing. You can provide a single method or a comma-separated list:
  - `DEG`: Differential Gene Expression analysis.
  - `WGCNA`: Weighted Gene Co-expression Network Analysis.
  - `DEG,WGCNA`: Runs both analysis workflows in parallel.

### Example Parameters

**Running both DEG and WGCNA with Kallisto:**

```bash
--pseudo_aligner kallisto --analysis_method DEG,WGCNA
```

**Running only WGCNA with HISAT2:**

```bash
--aligner hisat2 --analysis_method WGCNA
```

**Transcript-level DEG with Salmon:**

```bash
--pseudo_aligner salmon --analysis_method DEG
```

---

## Germline Variant Calling (GVC)

When using the HISAT2 genomic alignment route, the pipeline can perform high-sensitivity variant calling to identify SNPs and Indels. This workflow integrates the GATK4 engine and SnpEff for functional annotation.

### Automated Reference Management

The pipeline can automatically fetch and index necessary reference files based on your `--species` and `--genome` parameters.

- **SnpEff Annotation**: Automatically provided for all supported species/genome combinations.
- **dbSNP**: Automatically provided for most species **EXCEPT**:
  - `fruitfly` (`BDGP6`)
  - `worm` (`WBcel235`)
  - `monkey` (`Mmul_10`)
- **Indexing**: If you do not provide pre-built index files, the pipeline will automatically generate the `.fai`, `.dict`, and `interval_list` from your reference FASTA.

### Required & Optional Parameters

If you are working with a species listed in the exceptions above, or wish to use custom versions, provide the following:

- `--fai`: (Optional) Path to Samtools Fasta Index (`.fai`).
- `--dict`: (Optional) Path to GATK Sequence Dictionary (`.dict`).
- `--interval_list`: (Optional) Path to GATK Interval List (`.interval_list`).
- `--dbsnp`: Path to a dbSNP VCF file. (Required for GVC if not automatically supported).
- `--dbsnp_tbi`: (Optional) Path to the pre-built index (`.tbi`) for the dbSNP VCF.
- `--snpeff_genome`: The specific SnpEff database name (e.g., `GRCh38.105`).
- `--snpeff_db`: (Optional) Path to a pre-built SnpEff cache directory.
- `--knownsites`: Path to VCF files of known polymorphic sites (Required for BQSR).

### GVC Step-by-Step Parameters

#### 1. Interval Splitting & Parallelization

To optimize performance, the pipeline splits the genome into manageable segments.

- `--skip_interval_splitting`: Set to `true` to disable parallelization across genomic intervals.

#### 2. Duplicate Marking (Picard)

- `remove_duplicates`: Set to `true` to physically remove duplicate reads from the BAM (Default: `false`).
- `skip_picard_markduplicates`: Set to `true` to skip this step entirely.

#### 3. Base Quality Score Recalibration (BQSR)

- `skip_baserecalibration`: Set to `true` to skip.

> [!CAUTION]
> If you do not provide `--knownsites`, you **must** set `--skip_baserecalibration true` to prevent a crash.

#### 4. Variant Calling (HaplotypeCaller)

The pipeline identifies SNPs and small Indels by performing local de-novo assembly around active regions.

- `--gatk_hc_call_conf`: Minimum confidence threshold for calling a variant (Default: `20`).

#### 5. Hard Filtering (VariantFiltration)

Applies standard GATK filters to the raw VCF to reduce false positives.

- `--skip_variantfiltration`: Set to `true` to skip filtering.
- **Filter Thresholds**:
  - `gatk_vf_window_size`: Window size for filtering clusters (Default: `35`).
  - `gatk_vf_cluster_size`: Number of SNPs within window to trigger filter (Default: `3`).
  - `gatk_vf_fs_filter`: Fisher Strand (FS) threshold (Default: `30.0`).
  - `gatk_vf_qd_filter`: QualByDepth (QD) threshold (Default: `2.0`).

#### 6. Functional Annotation (SnpEff)

The final VCF is annotated to predict the biological impact of the variants (e.g., missense, nonsense, or synonymous mutations).

- `--snpeff_genome`: Specific SnpEff database name.
- `--snpeff_db`: (Optional) Path to a pre-built SnpEff cache directory.

**Standard GVC run with HISAT2 for a supported species (e.g., Human):**

```bash
--species human --genome GRCh38 \
--aligner hisat2 \
--dbsnp 'refs/dbsnp.vcf.gz' \
--knownsites 'refs/mills.vcf.gz,refs/1000G.vcf.gz' \
--snpeff_genome 'GRCh38.105'
```

**GVC for an unsupported species (e.g., Fruitfly) with manual BQSR skip:**

```bash
--species fruitfly --genome BDGP6 \
--aligner hisat2 \
--dbsnp 'refs/custom_dbsnp.vcf.gz' \
--skip_baserecalibration true
```

---

## Alternative Splicing (rMATS)

The pipeline uses [`rMATS`](http://rnaseq-mats.sourceforge.net/) to detect differential splicing events such as skipped exons (SE), retained introns (RI), and alternative 5'/3' splice sites.

### Activation Logic

- **Requirement 1**: `--aligner hisat2` must be selected.
- **Requirement 2**: `--analysis_method` must contain `AS` (e.g., `--analysis_method AS` or `--analysis_method AS,DEG`).

### rMATS Parameters

- `--statoff`: Set to `true` to skip statistical analysis and only perform read counting (Default: `false`).
- `--novelss`: Set to `true` to allow the detection of novel splice sites not present in the GTF (Default: `false`).
- `--allowclipping`: Set to `true` to allow soft-clipping of reads during analysis (Default: `false`).
- `--individualcounts`: Set to `true` to output counts for each individual replicate (Default: `false`).

### Filtering Significant Events

The pipeline includes an automated filtering step to extract high-confidence splicing events from the raw rMATS output.

#### Filtering Logic

The pipeline applies a rigorous `awk`-based filter to the Junction Count (JC) files using the following thresholds:

1. **Significance**: FDR (`--pvalue_threshold`) $\le 0.05$.
2. **Biological Effect**: $|\Delta \psi|$ (`--delta_psi`) $\ge 0.1$.
3. **Coverage**: Average read counts across replicates must be $\ge 10$ to ensure reliability.
4. **Inclusion Level**: Average $\psi$ values must be between $0.05$ and $0.95$ (filtering out constitutive splicing).

#### Filter Parameters

- `--event_types`: Comma-separated list of events to analyze (e.g., `SE,RI,MXE,A5SS,A3SS`). (Default: `SE`).
- `--delta_psi`: The minimum change in Percent Spliced In ($\psi$) required (Default: `0.1`).
- `--pvalue_threshold`: Used as the FDR (False Discovery Rate) cutoff (Default: `0.05`).

### Visualization (Sashimi Plots)

For every event that passes the filtering criteria, the pipeline automatically generates a **Sashimi Plot** using `rmats2sashimiplot`.

- **Output**: Detailed `.pdf` or `.png` plots showing the read coverage and junction spanning reads for each condition, allowing for visual validation of the splicing event.

### Example Parameters

**Analyzing Skipped Exons and Retained Introns with strict filtering:**

```bash
--aligner hisat2 --pseudo_aligner null \
--event_types 'SE,RI' \
--delta_psi 0.2 --pvalue_threshold 0.01
```

**Allowing novel splice site detection:**

```bash
--aligner hisat2 --pseudo_aligner null --novelss true
```

---

## Differential Exon Usage (DEXSeq)

The pipeline utilizes [`DEXSeq`](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) to identify differential exon usage. This method is more granular than standard gene-level analysis, allowing for the detection of specific isoform changes or alternative splicing effects that shift the expression of individual exons.

### Activation & Reference Requirements

- **Activation**: Triggered when `--analysis_method` contains `DEU` (e.g., `--analysis_method DEU,DEG`).
- **Requirement**: `--aligner hisat2` must be used.
- **GFF File**: DEXSeq requires a flattened GFF file.
  - The pipeline can automatically build this from your provided `--gtf`.
  - Alternatively, you can provide a pre-built file via `--gff`.

### DEXSeq Parameters

#### 1. Counting Phase

- `--aggregation`: Set to `true` to aggregate overlapping genes into a single gene cluster (Default: `false`).
- `--alignment_quality`: Minimum MAPQ score for a read to be included in the counting (Default: `20`).

#### 2. Statistical Analysis

- `--min_exonlength`: Exons shorter than this length (in bp) will be filtered out (Default: `18`).
- `--pvalue_threshold`: FDR (Adjusted p-value) cutoff for significance (Default: `0.05`).
- `--logfc_threshold`: Minimum absolute log2 fold change required for an exon to be considered significant (Default: `1`).

### Functional Enrichment (Optional)

If `DEU` is selected, you can optionally perform functional enrichment analysis on the genes associated with significant differential exon usage.

- `--enrichment_method`: Choose the enrichment type. Options: `GO`, `KEGG`, or `GO,KEGG`.
- **Logic**: The pipeline extracts the unique gene IDs from the significant DEU results and identifies enriched biological processes (GO) or metabolic pathways (KEGG).

### Example Parameters

**Running DEU with GO and KEGG enrichment:**

```bash
--aligner hisat2 --analysis_method DEU \
--gtf 'refs/genes.gtf' --enrichment_method GO,KEGG \
--pvalue_threshold 0.01 --logfc_threshold 1.5
```

**Combining DEU with standard DEG:**

```bash
--aligner hisat2 --analysis_method DEU,DEG \
--gtf 'refs/genes.gtf' --alignment_quality 30
```

---

## Differential Gene Expression (DEG)

The DEG module performs statistical testing to identify genes with significant expression changes across conditions or time points. Once identified, the pipeline can perform functional enrichment to provide biological context.

### Activation & Method Selection

- **Activation**: Include `DEG` in the `--analysis_method` parameter.
- **Method Choice**: Defined via `--diffexpr_method`.
  - `deseq2`: (**Default**) Optimized for small sample sizes; uses shrinkage estimation for dispersions.
  - `limma`: Best for large cohorts; utilizes linear modeling for high-speed computation.
  - `masigpro`: Specifically designed for **time-course** transcriptomic data.

### 1. Standard DEG Parameters (DESeq2 / Limma)

- `--pvalue_threshold`: FDR (adjusted p-value) cutoff for significance (Default: `0.05`).
- `--logfc_threshold`: Minimum absolute log2 fold change required (Default: `1`).
- `--ntop_genes`: The number of top-ranked significant genes to export into the **results CSV files** (Default: `100`).

### 2. Time-Course Parameters (maSigPro)

- **Samplesheet Requirement**: You must include a `timepoint` column as the **last column** in your `samplesheet.csv`. This column must be **numeric only** (e.g., `6`, `12`).
- `--rsq_threshold`: R-squared threshold to filter for genes that fit the temporal model (Default: `0.7`).
- `--cluster_method`: Method for grouping temporal profiles. Options: `hclust` (Default), `Mclust`, or `kmeans`.

### Workflow & Visualization

#### Standard Analysis (DESeq2 / Limma)

- **Quality Control**: PCA plots and Violin plots comparing **Non-filtered**, **Filtered**, and **Normalized** counts.
- **Statistical Plots**: Volcano plots and Heatmaps of the significant genes.
- **Limma Exclusive**: A Multi-Dimensional Scaling (MDS) plot is automatically generated.

#### Time-Course Analysis (maSigPro)

- **Cluster Heatmaps**: Genes are grouped into clusters based on their temporal patterns, with a dedicated heatmap generated for **each significant cluster**.
- **Quality Control**: Standard PCA and Violin plots (Raw/Filtered/Normalized).

### Functional Enrichment (Optional)

Triggered by the `--enrichment_method` parameter.

#### GO & KEGG Analysis

Available for all DEG methods. It identifies significant biological themes using enrichment p-values.

- `--ntop_processes`: Number of top-ranked pathways/terms to display in final plots (Default: `10`).
- **Directional Logic (DESeq2/Limma)**: Enrichment is performed separately for **Up-regulated** and **Down-regulated** significant genes.
- **Temporal Logic (maSigPro)**: Enrichment is performed independently for the significant genes within **each individual temporal cluster**.
- **GO Specifics**: Results are automatically split into **BP** (Biological Process), **MF** (Molecular Function), and **CC** (Cellular Component).

#### MSigDB GSEA (DESeq2 / Limma Only)

Uses a ranked list of all genes to identify broader biological signatures.

- `--msigdb_categories`: Comma-separated list of collections (e.g., `H,C2,C5`). (Default: `C1`).
- `--nes_threshold`: Minimum absolute Normalized Enrichment Score (Default: `1`).
- `--padj_gsea`: Adjusted p-value cutoff for GSEA results (Default: `0.25`).
- `--rank_method`: Method to order genes. Options: `logfc` (Default), `t_stat`, or `signed_significance` ($Log2FC \times -\log_{10}(p-value)$).

### Example Parameters

**Standard DEG with DESeq2 and strict GSEA:**

```bash
--analysis_method DEG --diffexpr_method deseq2 \
--enrichment_method MSigDB --msigdb_categories 'H,C5' \
--rank_method signed_significance --nes_threshold 1.5
```

**Time-course study using maSigPro with K-means clustering:**

```bash
--analysis_method DEG --diffexpr_method masigpro \
--cluster_method kmeans --rsq_threshold 0.8 --enrichment_method GO
```

---

## Weighted Gene Co-expression Network Analysis (WGCNA)

The WGCNA module identifies clusters of co-expressed genes and summarizes them using "module eigengenes." It helps identify key drivers (hub genes) associated with specific conditions or timepoints. Once identified, the pipeline can perform functional enrichment to provide biological context.

### Activation & Core Logic

- **Activation**: Include `WGCNA` in the `--analysis_method` parameter.
- **Time-Course Support**: Like `maSigPro`, WGCNA can correlate gene modules with the `timepoint` column provided in your samplesheet.

### WGCNA Network Parameters

- `--networktype`: Defines how the adjacency matrix is calculated. Options: `signed` (Default), `unsigned`, or `signed hybrid`.
- `--tomtype`: Type of Topological Overlap Matrix. Options: `signed` (Default), `none`, `signed 2`, `unsigned`, `unsigned 2`, `signed Nowick`, `signed Nowick 2`.
- `--sft_r2_threshold`: The R-squared threshold for picking the soft-thresholding power (Default: `0.8`).
- `--minmodulesize`: Minimum number of genes required to form a module (Default: `30`).
- `--deepsplit`: Sensitivity for module detection (Default: `2`).
- `--mergecutheight`: Dendrogram cut height for merging similar modules (Default: `0.25`).
- `--reassignthreshold`: p-value threshold for reassigning genes between modules (Default: `0`).

### Hub Gene & PPI Analysis

The pipeline identifies "Hub Genes" (the most connected genes in a module) based on Gene Significance (GS) and Module Membership (MM).

- `--min_gs`: Minimum Gene Significance to be considered a hub (Default: `0.5`).
- `--min_mm`: Minimum Module Membership (kME) to be considered a hub (Default: `0.8`).
- `--ntop_hubgenes`: The number of top-ranked hub genes to display in heatmaps and export to CSV (Default: `10`).
- `--cor_threshold`: Correlation threshold for network edge visualization (Default: `0.5`).

#### STRING PPI Integration

For significant modules, the pipeline uses the hub genes to query the **STRING database** for Protein-Protein Interactions (PPI).

- `--score_threshold`: Minimum interaction score for STRING edges (Default: `400`).

### Workflow & Visualization

- **Quality Control**: Standard PCA and Violin plots (Non-filtered vs. Filtered vs. Normalized).
- **Network Plots**: Cluster dendrograms (Module Colors) and Module-Trait heatmaps.
- **Hub Gene Visuals**: Heatmaps and CSV tables for the top hub genes in each module.
- **PPI Networks**: Visualization of protein interactions within significant modules.

### Functional Enrichment (Optional)

Triggered by the `--enrichment_method` parameter.

- **Logic**: Enrichment (GO/KEGG) is performed independently for the **significant hub genes** identified in each significant module.
- **Outputs**: Identifies the biological processes (BP, MF, CC) and pathways associated with specific network hubs.

### Example Parameters

**WGCNA with signed network and STRING PPI:**

```bash
--analysis_method WGCNA --networktype signed --tomtype signed \
--min_gs 0.6 --min_mm 0.85 --score_threshold 500
```

**WGCNA with GO/KEGG enrichment for hub genes:**

```bash
--analysis_method WGCNA --enrichment_method GO,KEGG \
--ntop_hubgenes 15 --mergecutheight 0.2
```

---

## Differential Isoform Usage (DIU) & Isoform Switching

The pipeline uses [`IsoformSwitchAnalyzeR`](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html) to identify functional shifts between transcript isoforms. This module moves beyond gene-level expression to find genes that "switch" their primary transcript between conditions. Once identified, the pipeline can perform functional enrichment to provide biological context.

### Activation & Essential Requirements

- **Activation**: Include `DIU` (for isoform switching) and/or `AS` (for alternative splicing analysis) in the `--analysis_method` parameter.
- **Requirement**: You **must** use a pseudo-aligner (`--pseudo_aligner salmon` or `kallisto`). This module is not compatible with the HISAT2 genomic alignment route.
- **Reference Files**:
  - `--transcriptome`: Reference cDNA sequences.
  - `--gtf_isoformswitchanalyzer`: A specific GTF file.
    - **Human & Mouse**: It is highly recommended to use **Ensembl** (including haplotype/scaffold annotations) or **GENCODE**.
    - **Other Species**: Standard Ensembl GTFs are typically sufficient.
    - **Note**: The Transcriptome and GTF versions **must** match exactly.

### Analysis Parameters

- `--pvalue_threshold`: FDR (adjusted p-value) cutoff for identifying significant isoform switches (Default: `0.05`).
- `--dif_cutoff`: The minimum absolute change in **dIF** (differential Isoform Fraction) required for a switch to be significant (Default: `0.1`).
- `--ntop_isoforms`: The number of top-ranked switching events to visualize in plots and export to the results CSV (Default: `10`).

### Workflow & Visualization

The pipeline identifies significant "switches" where a gene changes its primary isoform structure (e.g., from a long protein-coding transcript to a short non-coding one).

1. **Switch Detection**: Identifies genes where the relative contribution of isoforms significantly changes between conditions.
2. **Splicing Characterization**: If `AS` is selected, the pipeline identifies the specific splicing events (e.g., Exon Skipping, Alternative 3' Splice Sites) causing the switch.
3. **Isoform Switch Plots**: Generates detailed visualizations for the top `ntop_isoforms`, showing the transcript structure, domains, and expression shifts.

### Functional Enrichment (Optional)

Triggered by the `--enrichment_method` parameter.

- **Logic**: The pipeline extracts the parent genes of all significant isoform switching events and performs **GO** and **KEGG** enrichment.
- **Purpose**: Helps determine if isoform switching is impacting specific biological pathways or molecular functions.

### Example Parameters

**Running DIU with Salmon and GO enrichment:**

```bash
--pseudo_aligner salmon --analysis_method DIU --enrichment_method GO
```

**Analyzing both DIU and AS with strict thresholds:**

```bash
--pseudo_aligner kallisto --analysis_method DIU,AS \
--pvalue_threshold 0.01 --dif_cutoff 0.2 --ntop_isoforms 20
```

---

## Quality Control & Reporting

The pipeline performs automated quality assessment at every stage of the workflow. All metrics are aggregated into a single, interactive [`MultiQC`](https://multiqc.info/) report.

### QC Parameters

- `--skip_fastqc`: Set to `true` to skip **all** FastQC processes (Raw and Post-Trimming).
- `--skip_rseqc`: Set to `true` to skip all alignment-level quality control.
- `--rseqc_modules`: (Optional) A comma-separated list of specific modules to run.
- `--skip_multiqc`: Set to `true` to skip the final report aggregation.

> [!IMPORTANT]
> Some RSeQC modules, such as `genebody_coverage` and `tin`, can be computationally expensive and may increase memory usage for large datasets.

### QC Checkpoints

#### 1. Read Quality (FastQC)

Initial assessment of raw sequencing data and verification of adapter trimming efficiency using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

#### 2. Alignment Quality (RSeQC)

Post-alignment validation using [`RSeQC`](http://rseqc.sourceforge.net/). Supported modules include:

| Module              | Description                                                          |
| :------------------ | :------------------------------------------------------------------- | --- | --------------------- | ----------------------------------------------------- |
| `bam_stat`          | Summarizes mapping statistics and quality flags.                     |
| `genebody_coverage` | Checks for 5'/3' bias across gene bodies.                            |
| `infer_experiment`  | Verifies if the library is stranded (forward/reverse).               |
| `inner_distance`    | **(Paired-end only)** Calculates the insert size between read pairs. |     | `junction_annotation` | Compares detected splice junctions against reference. |
| `read_distribution` | Calculates % of reads in CDS, UTR, and Introns.                      |
| `read_duplication`  | Measures PCR and optical duplication levels.                         |
| `tin`               | **Transcript Integrity Number**: Assesses RNA degradation.           |

#### 3. Final Reporting (MultiQC)

The final report is generated in your output directory: `results/multiqc/multiqc_report.html`.

### Example Parameters

**Running only specific RSeQC modules:**

```bash
--rseqc_modules 'bam_stat,read_distribution,tin'
```

**Skipping all QC and Reporting for a fast run:**

```bash
--skip_fastqc --skip_multiqc
```

---

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run LoaiEletr/rnaseq \
   -profile docker \
   --input samplesheet.csv \
   --outdir ./results \
   --species human \
   --genome GRCh38 \
   --pseudo_aligner kallisto \
   --analysis_method DEG,WGCNA
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run LoaiEletr/rnaseq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
species: 'human'
genome: 'GRCh38'
pseudo_aligner: 'kallisto'
analysis_method: 'DEG,WGCNA'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull LoaiEletr/rnaseq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [LoaiEletr/rnaseq releases page](https://github.com/LoaiEletr/rnaseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
