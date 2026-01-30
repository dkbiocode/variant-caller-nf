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
  def samples_input = Channel
      .fromPath(params.samples_csv)
      .splitCsv(header: true)
      .map { row ->
          tuple(
              row.patient.trim(),
              row.tumor_srr.trim(),
              row.normal_srr.trim(),
              row.sample_type.trim()
          )
      }

  // Use multiMap to split into two channels without consuming the source
  def sample_channels = samples_input.multiMap { patient, tumor_srr, normal_srr, sample_type ->
      // Channel for pairing (preserve all metadata)
      pairing: tuple(patient, tumor_srr, normal_srr, sample_type)
      // Channel for extracting SRR IDs
      srr_extraction: tuple(tumor_srr, normal_srr)
  }

  // Channels created by multiMap
  samples_ch = sample_channels.pairing
  srr_pairs = sample_channels.srr_extraction

  // Extract unique SRR IDs for processing (both tumor and normal)
  all_srrs = srr_pairs
      .flatMap { tumor, normal -> [tumor, normal] }
      .unique()


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
    GATK(MARKDUPS.out.markdup_bam, DOWNLOAD_HG19.out.hg19.collect(), INDEX_FASTA.out.indexed_fasta.collect(), DOWNLOAD_GATK_RESOURCES.out.vcfs.collect(), DOWNLOAD_GATK_RESOURCES.out.indexes.collect())

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
    // IMPORTANT: Collect all final_bams first to avoid race condition
    // Convert to Map for easy lookup: SRR_ID -> [bam, bai]
    all_final_bams = final_bams
        .toList()
        .map { bam_list ->
            // Convert list of [srr, bam, bai] tuples to a Map
            def bam_map = [:]
            bam_list.each { tuple_item ->
                def srr = tuple_item[0]
                def bam = tuple_item[1]
                def bai = tuple_item[2]
                bam_map[srr] = [bam, bai]
            }
            return bam_map
        }

    paired_samples = samples_with_tumor_bam
        .combine(all_final_bams)  // Combine with the BAM map
        .map { item ->
            // Unpack: first 6 elements are from samples_with_tumor_bam, last is the bam map
            def normal_srr = item[0]
            def patient = item[1]
            def tumor_srr = item[2]
            def sample_type = item[3]
            def tumor_bam = item[4]
            def tumor_bai = item[5]
            def bam_map = item[6]  // This is a Map: SRR_ID -> [bam, bai]

            // Lookup the normal BAM from the map
            def normal_files = bam_map[normal_srr]

            tuple(
                patient,
                tumor_srr,
                normal_srr,
                sample_type,
                tumor_bam,
                tumor_bai,
                normal_files[0],  // normal_bam
                normal_files[1]   // normal_bai
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
        DOWNLOAD_HG19.out.hg19.collect(),
        INDEX_FASTA.out.indexed_fasta.collect()
    )

    VARSCAN_SOMATIC(SAMTOOLS_MPILEUP.out.pileups)

    VARSCAN_PROCESS(VARSCAN_SOMATIC.out.variants)

    // =========================================================================
    // SUMMARIZE AND AGGREGATE RESULTS
    // =========================================================================
    SUMMARIZE_VARIANTS(VARSCAN_PROCESS.out.somatic_hc)

    // Collect all variant summaries and create combined table
    AGGREGATE_RESULTS(SUMMARIZE_VARIANTS.out.summary.map { patient, sample_type, tsv -> tsv }.collect())
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
    def dict_name = hg19.toString().replaceAll(/\.fa$/, '.dict')
    """
    samtools faidx ${hg19}
    gatk CreateSequenceDictionary -R ${hg19} -O ${dict_name}
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
    """
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
    bwa mem -t ${task.cpus} -M -T 30 ${params.bwa_index_prefix} ${trimmed_r1} ${trimmed_r2} | \\
    samtools view -h -f 2 > ${srr}.sam
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
    def mem_gb = task.memory ? task.memory.toGiga() - 2 : 8
    def tmp_dir = task.executor == 'slurm' && task.scratch ? "\${SLURM_SCRATCH}" : "."
    """
    samtools addreplacerg -r ID:${srr} -r SM:${srr} -r PL:ILLUMINA -o ${srr}.rg.bam ${bam}

    picard -Xmx${mem_gb}g MarkDuplicates \\
        -I ${srr}.rg.bam \\
        -O ${srr}.markdup.bam \\
        -M ${srr}.marked_dup_metrics.txt \\
        --VERBOSITY DEBUG \\
        --TMP_DIR ${tmp_dir}
    """
}

process GATK {
    tag "${srr}"
    label "gatk"
    publishDir "${params.outdir}/final_bams", mode: 'symlink', pattern: "*.bam*"

    input:
    tuple val(srr), path(markdup_bam)
    path(ref_fasta)
    path(fasta_indices)
    path(gatk_vcfs)
    path(gatk_vcf_indices)

    output:
    tuple val(srr), path("${srr}.final.bam"), path("${srr}.final.bam.bai"), emit: final_bam

    script:
    def temp_dir = task.executor == 'slurm' ? "\${SLURM_SCRATCH}" : "."
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
    path(ref_fasta)
    path(fasta_indices)

    output:
    tuple val(patient), val(tumor_srr), val(normal_srr), val(sample_type),
          path("${patient}_${sample_type}.normal.pileup"),
          path("${patient}_${sample_type}.tumor.pileup"), emit: pileups

    script:
    def prefix = "${patient}_${sample_type}"
    """
    # Generate pileup for normal sample
    samtools mpileup \\
        -q 20 -Q 20 \\
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

process SUMMARIZE_VARIANTS {
    tag "patient${patient}_${sample_type}"
    publishDir "${params.varscan_results}/summaries", mode: 'copy'

    input:
    tuple val(patient), val(sample_type),
          path(snp_vcf), path(indel_vcf)

    output:
    tuple val(patient), val(sample_type), path("${patient}_${sample_type}_variants.tsv"), emit: summary

    script:
    """
    summarize_variants.py ${patient} ${sample_type} ${snp_vcf} ${indel_vcf} ${patient}_${sample_type}_variants.tsv
    """
}

process AGGREGATE_RESULTS {
    publishDir "${params.varscan_results}", mode: 'copy'

    input:
    path(variant_tsvs)

    output:
    path("all_variants_summary.tsv"), emit: combined_summary

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from pathlib import Path

    # Read all TSV files and combine
    all_dfs = []
    for tsv_file in Path('.').glob('*_variants.tsv'):
        try:
            df = pd.read_csv(tsv_file, sep='\\t')
            if len(df) > 0:
                all_dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {tsv_file}: {e}")

    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)

        # Sort by patient, sample type, chromosome, position
        combined['chr_num'] = combined['chr'].str.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')
        combined['chr_num'] = pd.to_numeric(combined['chr_num'], errors='coerce')
        combined = combined.sort_values(['patient', 'sample_type', 'chr_num', 'pos'])
        combined = combined.drop('chr_num', axis=1)

        combined.to_csv('all_variants_summary.tsv', sep='\\t', index=False)

        print(f"\\n=== Overall Summary ===")
        print(f"Total variants across all samples: {len(combined)}")
        print(f"Somatic variants: {len(combined[combined['somatic_status'] == 'Somatic'])}")
        print(f"\\nVariants per patient:")
        print(combined.groupby('patient').size())
        print(f"\\nVariants per sample type:")
        print(combined.groupby(['patient', 'sample_type']).size())
    else:
        # Create empty file if no variants found
        pd.DataFrame().to_csv('all_variants_summary.tsv', sep='\\t', index=False)
        print("No variants found in any sample")
    """
}
