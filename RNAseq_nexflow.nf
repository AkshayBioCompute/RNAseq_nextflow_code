#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.work_dir = '/home/akshay/Akshay/RNAseq/'
params.ref_genome = "${params.work_dir}reference_genome.fa"
params.annotation_file = "${params.work_dir}annotation.gtf"

// Set the working directory
workDir "${params.work_dir}"

// Channel for FASTQ files
FASTQ_FILES = Channel.fromFilePairs("${params.work_dir}*.fastq.gz", flat: true)

// Process to generate MD5 checksums
process generateMD5 {
    input:
    path fastq_file

    output:
    file 'md5_checksums.txt'

    script:
    """
    md5sum ${fastq_file} > md5_checksums.txt
    """
}

// Process to unzip FASTQ files
process unzipFastq {
    input:
    path fastq_file

    output:
    path fastq_file.replace('.gz', '')

    script:
    """
    gunzip ${fastq_file}
    """
}

// Process to run FastQC
process runFastQC {
    input:
    path fastq_file

    output:
    path "QC_reports/${fastq_file}.html"

    script:
    """
    fastqc -o QC_reports/ ${fastq_file}
    """
}

// Process to run MultiQC
process runMultiQC {
    input:
    path fastqc_report

    output:
    path 'QC_reports/multiqc_report.html'

    script:
    """
    multiqc -o QC_reports/ QC_reports/
    """
}

// Process to trim paired-end reads using Fastp
process trimReads {
    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple path("trimmed_${sample}_R1.fastq"), path("trimmed_${sample}_R2.fastq"),
          path("QC_reports/fastp_${sample}.html"),
          path("QC_reports/fastp_${sample}.json")

    script:
    """
    fastp -i ${r1} -I ${r2} -o trimmed_${sample}_R1.fastq -O trimmed_${sample}_R2.fastq \
    -h QC_reports/fastp_${sample}.html -j QC_reports/fastp_${sample}.json --thread 8
    """
}

// Process to build Hisat2 index
process buildHisat2Index {
    input:
    path ref_genome

    output:
    file 'index_prefix.1.ht2'

    script:
    """
    hisat2-build ${ref_genome} index_prefix
    """
}

// Process to align reads with Hisat2
process hisat2Align {
    input:
    tuple val(sample), path(r1), path(r2)

    output:
    path "${sample}.sam"

    script:
    """
    hisat2 -x index_prefix -1 ${r1} -2 ${r2} -S ${sample}.sam
    """
}

// Process to convert SAM to BAM, sort, and generate statistics
process convertSamToBam {
    input:
    path sam_file

    output:
    path "${sam_file.replace('.sam', '_sorted.bam')}", path("${sam_file.replace('.sam', '_alignment_statistics.txt')}")

    script:
    """
    samtools view -bS ${sam_file} | samtools sort -o ${sam_file.replace('.sam', '_sorted.bam')}
    samtools flagstat ${sam_file.replace('.sam', '_sorted.bam')} > ${sam_file.replace('.sam', '_alignment_statistics.txt')}
    """
}

// Process to run featureCounts
process runFeatureCounts {
    input:
    path sorted_bam_files

    output:
    path 'gene_counts.txt'

    script:
    """
    featureCounts -T 8 -a ${params.annotation_file} -o gene_counts.txt -g gene_id -t exon -s 1 -p ${sorted_bam_files}
    """
}

// Workflow definition
workflow {
    // Create output directories
    mkdir 'QC_reports'

    // Execute the processes
    md5_checksums = FASTQ_FILES | generateMD5
    unzipped_fastq = FASTQ_FILES | unzipFastq
    fastqc_reports = unzipped_fastq | runFastQC
    multiqc_report = fastqc_reports.collect() | runMultiQC
    trimmed_reads = FASTQ_FILES | trimReads
    index = buildHisat2Index(params.ref_genome)
    aligned_sams = trimmed_reads | hisat2Align
    sorted_bams = aligned_sams | convertSamToBam
    runFeatureCounts(sorted_bams)
}
