#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sralist = 'samples_of_interest.txt'
params.outdir = 'results'
params.fastp_reports = "${params.outdir}/fastp-reports"
params.trimmed_reads = "${params.outdir}/trimmed"
params.hg19 = "${params.outdir}/hg19"
params.bwa_index_prefix = "hg19"
params.gatk_resources = "${params.outdir}/gatk_resources"

// Initial data channel- the SRA accession list.
sra_id_ch = Channel
    .fromPath(params.sralist)
    .splitText()
    .map { it.trim() }
    .filter { it }  // removes empty lines

workflow {

    // download and index reference genome
    DOWNLOAD_HG19()
    BWA_INDEX(DOWNLOAD_HG19.out.hg19)
    INDEX_FASTA(DOWNLOAD_HG19.out.hg19)
    

    // Download GATK resources
    DOWNLOAD_GATK_RESOURCES()

    // main pipeline
    // Main pipeline
    sra_id_ch | FASTERQ_DUMP | FASTP \
        | BWA_MEM(BWA_INDEX.out.index.collect()) \
        | SAMTOOLS_SORT \
        | MARKDUPS \
        | GATK(INDEX_FASTA.out.indexed_fasta.collect(), DOWNLOAD_GATK_RESOURCES.out.vcfs.collect()) \
        | VARSCAN

}

process DOWNLOAD_GATK_RESOURCES {
    tag "gatk-resources"
    label "download"
    storeDir "${params.gatk_resources}"

    output:
    path("*.vcf"), emit: vcfs
    path("*.vcf.idx"), emit: indexes

    script:
    """
    # Download from GATK bundle - adjust URLs to your source
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg19/v0/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg19/v0/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg19/v0/dbsnp_138.hg19.vcf
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg19/v0/dbsnp_138.hg19.vcf.idx
    """
}

process INDEX_FASTA {
    tag "fasta-index"
    label "basic"
    storeDir "${params.hg19}"

    input:
    path(hg19)

    output:
    tuple path(hg19), path("${hg19}.fai"), path("*.dict"), emit: indexed_fasta

    script:
    """
    samtools faidx ${hg19}
    gatk CreateSequenceDictionary -R ${hg19}
    """
}

process FASTERQ_DUMP {
    tag "${srr}"
    label 'fasterq_dump' // defines task.cpus and task.memory
     
    input:
    val srr
    
    output:
    tuple val(srr), path("raw_reads/${srr}_1.fastq"), path("raw_reads/${srr}_2.fastq"), emit: reads
    
    script:
    """
    module load sra-toolkit/3.0.0
    
    fasterq-dump -f \\
        --threads ${task.cpus} \\
        --mem ${task.memory.toGiga()}GB \\
        -O raw_reads \\
        --temp \$SLURM_SCRATCH \\
        ${srr}
    """
}

process DOWNLOAD_HG19 {
    tag "hg19"
    label 'basic'
    storeDir "${params.hg19}", mode: 'symlink'

    output:
    path("hg19.fa"), emit: hg19

    script:
    """
    wget https://hgdownload.gi.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    curl -s https://hgdownload.gi.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt | grep hg19.fa.gz > hg19.fa.gz.md5
    md5sum -c hg19.fa.gz.md5
    gunzip hg19.fa.gz
    """
}

process FASTP {
    tag "${srr}"
    label 'fastp'
    publishDir "${params.fastp_reports}", mode: 'copy', pattern: "*.{json,html}"
    publishDir "${params.trimmed_reads}", mode: 'symlink', pattern: "*.trimmed.fastq"
    
    input:
    tuple val(srr), path(r1), path(r2)

    output:
    tuple val(srr), path("${srr}_1.trimmed.fastq"), path("${srr}_2.trimmed.fastq"), emit: trimmed
    path("*.json"), emit: json
    path("*.html"), emit: html

    script:
    """
    module load miniforge/24.11.3-0
    conda activate variant-calling

    fastp --thread ${task.cpus} \\
    -i ${r1} \\
    -I ${r2} \\
    -o ${srr}_1.trimmed.fastq \\
    -O ${srr}_2.trimmed.fastq \\
    --json ${srr}.json \\
    --html ${srr}.html
    """
}

process BWA_INDEX {
    tag "hg19-index"
    label "bwa_index"

    storeDir "${params.hg19}"

    input:
    path(hg19) // fasta file

    output:
    tuple path(hg19), path("${params.bwa_index_prefix}.*"), emit: index
    

    script:
    """
    bwa index -p ${params.bwa_index_prefix} ${hg19}
    """
}

process BWA_MEM {
    tag "${srr}"
    label "bwa_mem"

    input:
    tuple val(srr), path(trimmed_r1), path(trimmed_r2)
    path(hg19)

    output:
    tuple val(srr), path("${srr}.sam"), emit: aligned_sam

    script:
    """
    bwa mem -t ${task.cpus} ${params.hg19}/${params.bwa_index_prefix} ${trimmed_r1} ${trimmed_r2} > ${srr}.sam
    """
}

process SAMTOOLS_SORT {
    tag "${srr}"
    label "samsort"

    input:
    tuple val(srr), path(sam)

    output:
    tuple val(srr), path("${srr}.sorted.bam"), emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} -m 3840M -T \${SLURM_SCRATCH}/${srr}.tmp -O bam -o ${srr}.sorted.bam ${sam}
    """
}

process MARKDUPS {
    tag "${srr}"
    label "markdups"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(srr), path(bam)

    output:
    tuple val(srr), path("${srr}.markdup.bam"), emit: markdup_bam

    script:
    """
    samtools addreplacerg -r ID:${srr} -r SM:${srr} -r PL:ILLUMINA -o ${srr}.rg.bam ${bam}

    picard MarkDuplicates -I ${srr}.rg.bam -O ${srr}.markdup.bam -M ${srr}.marked_dup_metrics.txt --TMP_DIR \${SLURM_SCRATCH}
    """
}

process GATK {
    tag "${srr}"
    label "gatk"

    input:
    tuple val(srr), path(markdup_bam)
    path(fasta_files)   // dependency signal
    path(gatk_vcfs)     // dependency signal

    output:
    tuple val(srr), path("${srr}.final.bam"), path("${srr}.final.bam.bai"), emit: final_bam

    script:
    """
    # Base Quality Score Recalibration
    gatk BaseRecalibrator \\
        --java-options "-Djava.io.tmpdir=\${SLURM_SCRATCH}" \\
        -R ${params.hg19}/hg19.fa \\
        -I ${markdup_bam} \\
        --known-sites ${params.gatk_resources}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \\
        --known-sites ${params.gatk_resources}/dbsnp_138.hg19.vcf \\
        -O ${srr}.recal_data.table --tmp-dir \${SLURM_SCRATCH}

    # Apply the BQSR corrections
    gatk ApplyBQSR \\
        -R ${params.hg19}/hg19.fa \\
        -I ${markdup_bam} \\
        --bqsr-recal-file ${srr}.recal_data.table \\
        -O ${srr}.final.bam --tmp-dir \${SLURM_SCRATCH}

    # Index final BAM
    samtools index ${srr}.final.bam

"""
}
