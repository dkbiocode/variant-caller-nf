## Variant caller for tumor progression

Treatment resistant tumors develop as a result of positive selection in the pathways which a pharmaceutical agent targets. An example is the epidermal growth factor receptor gene EGFR accumulating a mutation that results in constitutive activation, while anti-EGFR treatment is meant to interrupt signaling in that pathway. 

### Approach

To track the progression of resistance to treatment, the authors in Lim et al. (2021)[^1] use circulating tumor DNA to detect changes in the patients' tumor genotype, using somatic controls and multiple time points of liquid biopsy. 

### Workflow 

[This nextflow pipeline](reports/flowchart.pdf) calls variants on Illumina reads, and has modifications to track mutations in tumor biopsies after Lim et al. (2021)[^1], 
which follows patients with colorectal cancer who become resistant to anti-EGFR treatment.

### Features

* Branching/converging workflow compares PBMC control with cfDNA in timecourse
* HPC configuration (SLURM) scales to 1000s of samples
* Checks SNPs for select loci (targetted sequencing)
* Compares varscan and lofreq variant callers

### Results

* Biopsy resistance mutations confirmed
* Acquired resistance mutations tentative - needs calibration
* KRAS, MAP2K1, TP53 searched

[^1]: Circulating tumor DNA sequencing in colorectal cancer patients treated with first-line chemotherapy with anti-EGFR. [Scientific Reports volume 11, Article number: 16333 (2021) ](https://www.nature.com/articles/s41598-021-95345-4)
