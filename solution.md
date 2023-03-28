This pipeline can be used for the analysis of RNA-seq data:

- Downloads the reference genome and related files (Download_ref).
- Generates Bowtie2 indices for the reference genome (Bowtie_index).
- Trims adapter sequences and low-quality bases from the input reads (TrimGalore).
- Aligns the trimmed reads to the reference genome (Bowtie_Samtools).
- Counts the number of reads aligned to each gene in the reference genome (FeatureCounts).
- Adds gene names to the counts file (Add_gene_names).
- Generates quality control reports for the trimmed reads (MultiQC).

```mermaid
graph LR
    A[Download_ref] --> B[Bowtie_index]
    C{input_reads} --> D[TrimGalore]
    D --> E[MultiQC]
    B --> F[Bowtie_Samtools]
    D --> F
    F --> G[FeatureCounts]
    B --> G
    G --> H[Add_gene_names]
 ```
