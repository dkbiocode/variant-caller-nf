#!/usr/bin/env bash
# Quick script to check variant calling results on Alpine

echo "==================================================================="
echo "Checking Somatic Variant Calls for C_fw019"
echo "==================================================================="
echo ""

echo "SNP variants (excluding headers):"
for vcf in results/varscan/processed/*Somatic.hc.vcf; do
    if [[ -f "$vcf" ]]; then
        sample=$(basename "$vcf" | sed 's/.snp.Somatic.hc.vcf//' | sed 's/.indel.Somatic.hc.vcf//')
        count=$(grep -v "^#" "$vcf" | wc -l)
        echo "  $sample: $count variants"
    fi
done

echo ""
echo "==================================================================="
echo "Sample High-VAF Somatic SNPs (VAF > 5%):"
echo "==================================================================="
echo ""

for vcf in results/varscan/processed/*.snp.Somatic.hc.vcf; do
    if [[ -f "$vcf" ]]; then
        sample=$(basename "$vcf" .snp.Somatic.hc.vcf)
        echo "--- $sample ---"

        # Extract variants with FREQ > 5% from tumor sample
        grep -v "^#" "$vcf" | awk -F'\t' '{
            # Get FORMAT field indices
            split($9, format, ":");
            for (i in format) {
                if (format[i] == "FREQ") freq_idx = i;
            }

            # Get tumor FREQ value (11th column)
            split($11, tumor, ":");
            freq = tumor[freq_idx];
            gsub("%", "", freq);

            if (freq+0 > 5) {
                print $1 "\t" $2 "\t" $4 "\t" $5 "\t" freq "%";
            }
        }' | head -20
        echo ""
    fi
done

echo "==================================================================="
echo "Quick summary by sample type:"
echo "==================================================================="
echo ""

# Count by sample type
for type in biopsy plasma_baseline plasma_tp2 plasma_tp5; do
    snp_count=$(grep -v "^#" results/varscan/processed/C_fw019_${type}.snp.Somatic.hc.vcf 2>/dev/null | wc -l)
    indel_count=$(grep -v "^#" results/varscan/processed/C_fw019_${type}.indel.Somatic.hc.vcf 2>/dev/null | wc -l)
    echo "$type:"
    echo "  SNPs: $snp_count"
    echo "  INDELs: $indel_count"
    echo "  Total: $((snp_count + indel_count))"
    echo ""
done
