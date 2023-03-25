#!/usr/bin/env nextflow
/*
############################################

Your task if to draw a flow chart of the processes of this pipeline
You should also figure out what the individual prosesses do and what the pipeline can be used for. 

The 

############################################
*/
nextflow.enable.dsl=2
workflow {
    // Collect all read and trim them
    read_ch = Channel.fromPath(params.input_reads)
    .map { file -> tuple(file.baseName, file) }
    bowtie_index_ch = Bowtie_index(download_ref_ch)
    trim_ch = TrimGalore(read_ch) 
    MultiQC(trim_ch[1].collect().flatten()) 
    download_ref_ch = Download_ref()
    bowtie_samtools_ch = Bowtie_Samtools(bowtie_index_ch, trim_ch)
    featurecounts_ch = FeatureCounts(bowtie_samtools_ch[0].collect().flatten(), download_ref_ch.collect().flatten())
    Add_gene_names(featurecounts_ch.collect(), download_ref_ch.collect().flatten())

}

process MultiQC {

	publishDir "${params.directory_out}/data/RNA/MultiQC/"

	input:
	file reports

	output:
	file "multiqc_report.html"
	file "multiqc_data"

    script:
    """
    multiqc -f ${params.directory_out}/data/RNA/trimGalore/*
    """
}

process Download_ref {

    publishDir "${params.directory_out}/data/RNA/references/"  

    output:
    path "ncbi_dataset"

    script:
    """
    # Download referenses
    datasets download genome accession GCF_000005845.2,GCF_000750555.1 --include gff3,rna,cds,protein,genome,gtf,seq-report,gbff 
    # unzio downloaded file
    unzip ncbi_dataset.zip 
    """
}

process TrimGalore {

    publishDir "${params.directory_out}/data/RNA/trimGalore/"  
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*R1_001_trimmed.fq.gz")
    path "*fastqc*"
    path "*trimming_report.txt"


    script:
    """
    trim_galore --fastqc --adapter "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" --phred33 -e 0.1 --quality 28 --cores ${task.cpus} ${reads}
    """
}

process Bowtie_Samtools {

	publishDir "${params.directory_out}/data/RNA/Bowtie_samtools/"

	input:
	path reference_genome
    tuple val(sample_id), path(read)
    val x
    val y

	output:
	path "*.sorted.bam"
    path "*_.txt"

    script:
    """
     python ${params.scripts}/Bowtie_Samtools.py  --sample_id ${sample_id} --read ${read}
    """
}

process FeatureCounts {

	publishDir "${params.directory_out}/data/RNA/FeatureCounts/"

	input:
    path mapping
	path reference_genome
    
	output:

    path '*txt'

    script:

    """
    featureCounts -a ncbi_dataset/data/GCF_000005845.2/genomic.gtf \
    -t gene, \
    -g gene_id \
    -o featurecounts_comb.txt \
    -T ${task.cpus} \
     ${params.directory_out}/data/RNA/Bowtie_samtools/*sorted.bam 2> featurecounts_mapp_rate.txt
    """
}

process Add_gene_names {

    publishDir "${params.directory_out}/data/RNA/FeatureCounts/"  
    
    input:
    path file
    path references

    output:
    path "featurecounts_comb.txt"


    script:
    """
    python ${params.scripts}/Add_gene_name.py 
    """
}

process Bowtie_index {

	publishDir "${params.directory_out}/data/RNA/bowtie2_index/"

	input:
	file reports

	output:
	file "MG1655*"

    script:
    """
    bowtie2-build ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna MG1655
    """
}