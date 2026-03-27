## Variant caller

This nextflow pipeline calls variants on illumina reads, and has modifications to track mutations in tumor biopsies after Lim et al. (2021)[^1], 
which follows patients with colorectal cancer who become resistant to anti-EGFR treatment.

### Features

* Branching/converging workflow compares PBMC control with cfDNA in timecourse
* HPC configuration (SLURM) scales to 1000s of samples
* Checks SNPs for select loci (targetted sequencing)

### Results

* Biopsy resistance mutations confirmed
* Acquired resistance mutations tentative
* KRAS, MAP2K1, TP53 searched

[^1]: Circulating tumor DNA sequencing in colorectal cancer patients treated with first-line chemotherapy with anti-EGFR. [Scientific Reports volume 11, Article number: 16333 (2021) ](https://www.nature.com/articles/s41598-021-95345-4)
