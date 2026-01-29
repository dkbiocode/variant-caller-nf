#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.samples_csv = 'samples.csv'
params.outdir = 'results'
params.fastp_reports = "${params.outdir}/fastp-reports"
params.trimmed_reads = "${params.outdir}/trimmed"
params.hg19 = "${params.outdir}/hg19"
params.bwa_index_prefix = "hg19"
params.gatk_resources = "${params.outdir}/gatk_resources"
params.varscan_results = "${params.outdir}/varscan"

// VarScan parameters
params.min_coverage = 8
params.min_coverage_normal = 10
params.min_coverage_tumor = 6
params.min_var_freq = 0.10
params.min_freq_for_hom = 0.75
params.somatic_p_value = 0.05
params.strand_filter = 1

// Parse samples CSV to create tumor-normal pairs
// Each row contains: patient, tumor_srr, normal_srr, sample_type
samples_ch = Channel
    .fromPath(params.samples_csv)
    .splitCsv(header: true)
    .map { row -> 
        tuple(
            row.patient,
            row.tumor_srr,
            row.normal_srr,
            row.sample_type
        )
    }

// Extract unique SRR IDs for processing (both tumor and normal)
// We need to process each SRR only once even if it appears multiple times
tumor_srrs = samples_ch.map { patient, tumor, normal, type -> tumor }
normal_srrs = samples_ch.map { patient, tumor, normal, type -> normal }
all_srrs = tumor_srrs.mix(normal_srrs).unique()

workflow {

    // =========================================================================
    // REFERENCE GENOME AND RESOURCES
    // =========================================================================
    DOWNLOAD_HG19()
    BWA_INDEX(DOWNLOAD_HG19.out.hg19)
    INDEX_FASTA(DOWNLOAD_HG19.out.hg19)
    DOWNLOAD_GATK_RESOURCES()

    // =========================================================================
    // PREPROCESS ALL UNIQUE SAMPLES
    // =========================================================================
    // Process each unique SRR through: download -> trim -> align -> sort -> markdup -> BQSR
    FASTERQ_DUMP(all_srrs)
    FASTP(FASTERQ_DUMP.out.reads)
    BWA_MEM(FASTP.out.trimmed, BWA_INDEX.out.index.collect())
    SAMTOOLS_SORT(BWA_MEM.out.aligned_sam)
    MARKDUPS(SAMTOOLS_SORT.out.sorted_bam)
    GATK(MARKDUPS.out.markdup_bam, INDEX_FASTA.out.indexed_fasta.collect(), DOWNLOAD_GATK_RESOURCES.out.vcfs.collect())

    final_bam = GATK.out.final_bam

    // =========================================================================
    // PAIR TUMOR AND NORMAL SAMPLES BY PATIENT
    // =========================================================================
    // 
    // Goal: For each row in samples.csv, match the processed tumor BAM 
    //       with its corresponding normal BAM from the same patient.
    //
    // Input (samples_ch from CSV):
    //   [patient=7,  tumor_srr=SRR13974111, normal_srr=SRR13973950, type=biopsy]
    //   [patient=7,  tumor_srr=SRR13973861, normal_srr=SRR13973950, type=plasma]
    //   ...
    //
    // GATK output (final BAMs keyed by SRR):
    //   [SRR13974111, bam, bai]
    //   [SRR13973950, bam, bai]
    //   ...
    //
    // Desired output (paired_samples):
    //   [patient, tumor_srr, normal_srr, type, tumor_bam, tumor_bai, normal_bam, normal_bai]

    // Step 1: Create channel of final BAMs keyed by SRR ID
    // tuple(srr, bam, bai)
    final_bams = final_bam

    // Step 2: Join tumor BAMs to sample metadata
    // - Rekey samples_ch by tumor_srr for joining
    // - Join with final_bams to attach the tumor BAM files
    samples_with_tumor_bam = samples_ch
        .map { patient, tumor_srr, normal_srr, sample_type ->
            tuple(tumor_srr, patient, normal_srr, sample_type)
        }
        .join(final_bams)  // joins on first element (tumor_srr)
        .map { tumor_srr, patient, normal_srr, sample_type, tumor_bam, tumor_bai ->
            tuple(
                normal_srr,   // rekey by normal_srr for next join
                patient,
                tumor_srr,
                sample_type,
                tumor_bam,
                tumor_bai
            )
        }
    // Result: tuple(normal_srr, patient, tumor_srr, sample_type, tumor_bam, tumor_bai)

    // Step 3: Join normal BAMs to complete the pairing
    // - Join on normal_srr to attach normal BAM files
    // - Restructure to final format
    paired_samples = samples_with_tumor_bam
        .join(final_bams)  // joins on first element (normal_srr)
        .map { normal_srr, patient, tumor_srr, sample_type, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(
                patient,
                tumor_srr,
                normal_srr,
                sample_type,
                tumor_bam,
                tumor_bai,
                normal_bam,
                normal_bai
            )
        }
    // Result: tuple(patient, tumor_srr, normal_srr, sample_type, tumor_bam, tumor_bai, normal_bam, normal_bai)

    // Debug: uncomment to verify pairing
    // paired_samples.view { "PAIRED: patient=${it[0]} tumor=${it[1]} normal=${it[2]} type=${it[3]}" }

    // =========================================================================
    // SOMATIC VARIANT CALLING
    // =========================================================================
    SAMTOOLS_MPILEUP(
        paired_samples,
        INDEX_FASTA.out.indexed_fasta.collect()
    )

    VARSCAN_SOMATIC(SAMTOOLS_MPILEUP.out.pileups)

    VARSCAN_PROCESS(VARSCAN_SOMATIC.out.variants)
}

// =============================================================================
// REFERENCE AND RESOURCE DOWNLOAD PROCESSES
// =============================================================================

process DOWNLOAD_HG19 {
    tag "hg19"
    label 'basic'
    storeDir "${params.hg19}"

    output:
    path("hg19.fa"), emit: hg19

    script:
    """
    curl -O https://hgdownload.gi.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    curl -s https://hgdownload.gi.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt | grep hg19.fa.gz > hg19.fa.gz.md5
    md5sum -c hg19.fa.gz.md5
    gunzip hg19.fa.gz
    """
}

process DOWNLOAD_GATK_RESOURCES {
    tag "gatk-resources"
    label "download"
    storeDir "${params.gatk_resources}"

    output:
    path("*.vcf.gz"), emit: vcfs
    path("*.vcf.gz.tbi"), emit: indexes

    script:
    """
    # Download Mills indels from Broad FTP (requires special login)
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

    # Re-compress with bgzip and index (file is gzipped but not bgzipped)
    gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
    bgzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
    tabix -p vcf Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

    # Download dbSNP from NCBI latest release (GCF_000001405.25 = GRCh37/hg19)
    # Using latest_release which points to most recent build
    wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz -O dbsnp_138.hg19.vcf.gz
    wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi -O dbsnp_138.hg19.vcf.gz.tbi
    """
}

// =============================================================================
// INDEXING PROCESSES
// =============================================================================

process INDEX_FASTA {
    tag "fasta-index"
    label "basic"
    storeDir "${params.hg19}"

    input:
    path(hg19)

    output:
    tuple path("${hg19}.fai"), path("*.dict"), emit: indexed_fasta

    script:
    """
    samtools faidx ${hg19}
    gatk CreateSequenceDictionary -R ${hg19}
    """
}

process BWA_INDEX {
    tag "hg19-index"
    label "bwa_index"
    storeDir "${params.hg19}"

    input:
    path(hg19)

    output:
    path("hg19.*"), emit: index

    script:
    """
    bwa index -p ${params.bwa_index_prefix} ${hg19}
    """
}

// =============================================================================
// PREPROCESSING PROCESSES
// =============================================================================

process FASTERQ_DUMP {
    tag "${srr}"
    label 'fasterq_dump'
    storeDir "${params.outdir}/sra_downloads/${srr}"

    input:
    val srr

    output:
    tuple val(srr), path("${srr}_1.fastq"), path("${srr}_2.fastq"), emit: reads

    script:
    def temp_dir = task.executor == 'slurm' ? "\$SLURM_SCRATCH" : "\$TMPDIR"
    def module_load = task.executor == 'slurm' ? "module load sra-toolkit/3.0.0" : ""
    """
    ${module_load}

    fasterq-dump -f \\
        --threads ${task.cpus} \\
        --mem ${task.memory.toGiga()}GB \\
        --temp ${temp_dir} \\
        ${srr}
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
    def conda_init = task.executor == 'slurm' ? """
    module load miniforge
    conda activate variant-calling
    """ : ""
    """
    ${conda_init}

    fastp --thread ${task.cpus} \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${srr}_1.trimmed.fastq \\
        -O ${srr}_2.trimmed.fastq \\
        --json ${srr}.json \\
        --html ${srr}.html
    """
}

process BWA_MEM {
    tag "${srr}"
    label "bwa_mem"

    input:
    tuple val(srr), path(trimmed_r1), path(trimmed_r2)
    path(index_files)

    output:
    tuple val(srr), path("${srr}.sam"), emit: aligned_sam

    script:
    """
    bwa mem -t ${task.cpus} ${params.bwa_index_prefix} ${trimmed_r1} ${trimmed_r2} > ${srr}.sam
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
    publishDir "${params.outdir}/markdup_metrics", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(srr), path(bam)

    output:
    tuple val(srr), path("${srr}.markdup.bam"), emit: markdup_bam
    path("*.txt"), emit: metrics

    script:
    """
    module load miniforge
    conda activate variant-calling


    samtools addreplacerg -r ID:${srr} -r SM:${srr} -r PL:ILLUMINA -o ${srr}.rg.bam ${bam}

    picard MarkDuplicates \\
        -I ${srr}.rg.bam \\
        -O ${srr}.markdup.bam \\
        -M ${srr}.marked_dup_metrics.txt \\
        --TMP_DIR \${SLURM_SCRATCH}
    """
}

process GATK {
    tag "${srr}"
    label "gatk"
    publishDir "${params.outdir}/final_bams", mode: 'symlink', pattern: "*.bam*"

    input:
    tuple val(srr), path(markdup_bam)
    path(fasta_files)
    path(gatk_vcfs)

    output:
    tuple val(srr), path("${srr}.final.bam"), path("${srr}.final.bam.bai"), emit: final_bam

    script:
    def temp_dir = task.executor == 'slurm' ? "\${SLURM_SCRATCH}" : "."
    // Find the reference fasta - it should be named hg19.fa
    def ref_fasta = "hg19.fa"
    // Find the known sites VCFs
    def mills_vcf = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
    def dbsnp_vcf = "dbsnp_138.hg19.vcf.gz"
    """
    # Base Quality Score Recalibration
    gatk BaseRecalibrator \\
        --java-options "-Djava.io.tmpdir=${temp_dir}" \\
        -R ${ref_fasta} \\
        -I ${markdup_bam} \\
        --known-sites ${mills_vcf} \\
        --known-sites ${dbsnp_vcf} \\
        -O ${srr}.recal_data.table \\
        --tmp-dir ${temp_dir}

    # Apply the BQSR corrections
    gatk ApplyBQSR \\
        -R ${ref_fasta} \\
        -I ${markdup_bam} \\
        --bqsr-recal-file ${srr}.recal_data.table \\
        -O ${srr}.final.bam \\
        --tmp-dir ${temp_dir}

    # Index final BAM
    samtools index ${srr}.final.bam
    """
}

// =============================================================================
// SOMATIC VARIANT CALLING PROCESSES
// =============================================================================

process SAMTOOLS_MPILEUP {
    tag "patient${patient}_${sample_type}"
    label "mpileup"

    input:
    tuple val(patient), val(tumor_srr), val(normal_srr), val(sample_type), 
          path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(fasta_files)

    output:
    tuple val(patient), val(tumor_srr), val(normal_srr), val(sample_type),
          path("${patient}_${sample_type}.normal.pileup"), 
          path("${patient}_${sample_type}.tumor.pileup"), emit: pileups

    script:
    def prefix = "${patient}_${sample_type}"
    def ref_fasta = "hg19.fa"
    """
    # Generate pileup for normal sample
    samtools mpileup \\
        -q 1 \\
        -f ${ref_fasta} \\
        ${normal_bam} \\
        > ${prefix}.normal.pileup

    # Generate pileup for tumor sample
    samtools mpileup \\
        -q 1 \\
        -f ${ref_fasta} \\
        ${tumor_bam} \\
        > ${prefix}.tumor.pileup
    """
}

process VARSCAN_SOMATIC {
    tag "patient${patient}_${sample_type}"
    label "varscan"
    publishDir "${params.varscan_results}/raw", mode: 'copy'

    input:
    tuple val(patient), val(tumor_srr), val(normal_srr), val(sample_type),
          path(normal_pileup), path(tumor_pileup)

    output:
    tuple val(patient), val(tumor_srr), val(normal_srr), val(sample_type),
          path("${patient}_${sample_type}.snp.vcf"),
          path("${patient}_${sample_type}.indel.vcf"), emit: variants

    script:
    def prefix = "${patient}_${sample_type}"
    """
    module load miniforge
    conda activate variant-calling

    varscan somatic \\
        ${normal_pileup} \\
        ${tumor_pileup} \\
        ${prefix} \\
        --min-coverage ${params.min_coverage} \\
        --min-coverage-normal ${params.min_coverage_normal} \\
        --min-coverage-tumor ${params.min_coverage_tumor} \\
        --min-var-freq ${params.min_var_freq} \\
        --min-freq-for-hom ${params.min_freq_for_hom} \\
        --somatic-p-value ${params.somatic_p_value} \\
        --strand-filter ${params.strand_filter} \\
        --output-vcf 1
    """
}

process VARSCAN_PROCESS {
    tag "patient${patient}_${sample_type}"
    label "varscan"
    publishDir "${params.varscan_results}/processed", mode: 'copy'

    input:
    tuple val(patient), val(tumor_srr), val(normal_srr), val(sample_type),
          path(snp_vcf), path(indel_vcf)

    output:
    tuple val(patient), val(sample_type),
          path("${patient}_${sample_type}.snp.Somatic.hc.vcf"),
          path("${patient}_${sample_type}.indel.Somatic.hc.vcf"), emit: somatic_hc
    tuple val(patient), val(sample_type),
          path("${patient}_${sample_type}.snp.Germline.hc.vcf"),
          path("${patient}_${sample_type}.indel.Germline.hc.vcf"), emit: germline_hc
    tuple val(patient), val(sample_type),
          path("${patient}_${sample_type}.snp.LOH.hc.vcf"),
          path("${patient}_${sample_type}.indel.LOH.hc.vcf"), emit: loh_hc

    script:
    def prefix = "${patient}_${sample_type}"
    """
    module load miniforge
    conda activate variant-calling

    # Process SNPs - isolate calls by type and confidence
    varscan processSomatic ${snp_vcf} \\
        --min-tumor-freq ${params.min_var_freq} \\
        --max-normal-freq 0.05 \\
        --p-value ${params.somatic_p_value}

    # Process Indels - isolate calls by type and confidence
    varscan processSomatic ${indel_vcf} \\
        --min-tumor-freq ${params.min_var_freq} \\
        --max-normal-freq 0.05 \\
        --p-value ${params.somatic_p_value}
    """
}
