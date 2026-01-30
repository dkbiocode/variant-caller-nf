  # Comprehensive analysis
  echo "=== Total variants ==="
  wc -l results/varscan/all_variants_summary.tsv

  echo -e "\n=== Variants per patient ==="
  tail -n +2 results/varscan/all_variants_summary.tsv | cut -f1 | sort | uniq -c

  echo -e "\n=== Somatic vs other ==="
  tail -n +2 results/varscan/all_variants_summary.tsv | cut -f8 | sort | uniq -c

  echo -e "\n=== High-VAF somatic variants (>10%) ==="
  awk -F'\t' '$8 == "Somatic" && $7 > 10 {print $1, $2":"$3, $4">"$5, "VAF="$7"%"}' results/varscan/all_variants_summary.tsv | head -20

  echo -e "\n=== Checking for resistance mutations ==="
  ./check_resistance_mutations.sh
