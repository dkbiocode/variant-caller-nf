#!/usr/bin/env bash
set -e

cd results/lofreq/summaries 
# Grab header from first file, then append all files without their headers
head -1 $(ls *_lofreq_variants.tsv | head -1) > ../all_lofreq_variants.tsv
tail -n +2 -q *_lofreq_variants.tsv >> ../all_lofreq_variants.tsv
