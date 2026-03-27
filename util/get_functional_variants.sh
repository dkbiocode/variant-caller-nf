#!/usr/bin/env bash -ue

echo '# Functional variants only (missense, nonsense, frameshift, etc.)'
awk -F'\t' '$8=="KRAS" && ($9 ~ /missense|nonsense|frameshift|stop_gained|stop_lost|start_lost/)' results/varscan/all_variants_summary.tsv

awk -F'\t' '$8=="MAP2K1" && ($9 ~ /missense|nonsense|frameshift|stop_gained|stop_lost|start_lost/)' results/varscan/all_variants_summary.tsv

awk -F'\t' '$8=="TP53" && ($9 ~ /missense|nonsense|frameshift|stop_gained|stop_lost|start_lost/)' results/varscan/all_variants_summary.tsv

echo '2. What ARE all these variants?'

echo '# Breakdown by variant effect type'
echo "=== KRAS variant types ==="
awk -F'\t' '$8=="KRAS" {print $9}' results/varscan/all_variants_summary.tsv | sort | uniq -c | sort -rn

echo "=== MAP2K1 variant types ==="
awk -F'\t' '$8=="MAP2K1" {print $9}' results/varscan/all_variants_summary.tsv | sort | uniq -c | sort -rn

echo "=== TP53 variant types ==="
awk -F'\t' '$8=="TP53" {print $9}' results/varscan/all_variants_summary.tsv | sort | uniq -c | sort -rn

echo '3. Which are HIGH impact?'

# High-impact variants in resistance genes
awk -F'\t' '($8=="KRAS" || $8=="MAP2K1" || $8=="TP53") && $10=="HIGH"' results/varscan/all_variants_summary.tsv


