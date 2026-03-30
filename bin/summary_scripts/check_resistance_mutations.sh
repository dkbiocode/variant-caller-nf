#!/usr/bin/env bash
# Check for known resistance mutations from Table 4 of the publication

echo "==================================================================="
echo "Searching for Known Resistance Mutations (from Table 4)"
echo "==================================================================="
echo ""

# Define known resistance positions (hg19 coordinates)
declare -A RESISTANCE_SITES=(
    ["chr12:25398284"]="KRAS c.35G>T p.Gly12Val"
    ["chr12:25398285"]="KRAS c.34G>T p.Gly12Cys"
    ["chr12:25380275"]="KRAS c.183A>C p.Gln61His"
    ["chr15:66727455"]="MAP2K1 c.171G>T p.Lys57Asn"
    ["chr15:66727454"]="MAP2K1 c.170A>C p.Lys57Thr"
    ["chr17:7577548"]="TP53 c.733G>A p.Gly245Ser"
    ["chr5:112175639"]="APC c.4348C>T p.Arg1450*"
    ["chr5:112177901"]="APC c.6610C>T p.Arg2204*"
    ["chr10:89711900"]="PTEN c.518G>A p.Arg173His"
)

# Also check without "chr" prefix (VarScan might use different formats)
declare -A RESISTANCE_SITES_NOCHROM=(
    ["12:25398284"]="KRAS c.35G>T p.Gly12Val"
    ["12:25398285"]="KRAS c.34G>T p.Gly12Cys"
    ["12:25380275"]="KRAS c.183A>C p.Gln61His"
    ["15:66727455"]="MAP2K1 c.171G>T p.Lys57Asn"
    ["15:66727454"]="MAP2K1 c.170A>C p.Lys57Thr"
    ["17:7577548"]="TP53 c.733G>A p.Gly245Ser"
    ["5:112175639"]="APC c.4348C>T p.Arg1450*"
    ["5:112177901"]="APC c.6610C>T p.Arg2204*"
    ["10:89711900"]="PTEN c.518G>A p.Arg173His"
)

found_any=0

for vcf in results/varscan/processed/*.Somatic.hc.vcf; do
    if [[ ! -f "$vcf" ]]; then
        continue
    fi

    sample=$(basename "$vcf" | sed 's/.snp.Somatic.hc.vcf//' | sed 's/.indel.Somatic.hc.vcf//')

    # Check each known position
    for pos in "${!RESISTANCE_SITES[@]}"; do
        chr=$(echo "$pos" | cut -d: -f1)
        coord=$(echo "$pos" | cut -d: -f2)

        # Try with chr prefix
        match=$(grep -v "^#" "$vcf" | grep "^${chr}[[:space:]]${coord}[[:space:]]")

        # Try without chr prefix if not found
        if [[ -z "$match" ]]; then
            chr_num=$(echo "$chr" | sed 's/chr//')
            match=$(grep -v "^#" "$vcf" | grep "^${chr_num}[[:space:]]${coord}[[:space:]]")
        fi

        if [[ -n "$match" ]]; then
            echo "*** RESISTANCE MUTATION FOUND ***"
            echo "Sample: $sample"
            echo "Position: $pos"
            echo "Mutation: ${RESISTANCE_SITES[$pos]}"
            echo "VCF line:"
            echo "$match"
            echo ""
            found_any=1
        fi
    done
done

if [[ $found_any -eq 0 ]]; then
    echo "No exact matches to Table 4 resistance mutations found."
    echo ""
    echo "This could mean:"
    echo "  1. Patient C_fw019 doesn't have these specific mutations"
    echo "  2. Different chromosome naming (checking both formats)"
    echo "  3. Low coverage at these positions"
    echo "  4. VAF below detection threshold"
    echo ""
    echo "Let's check what's in the KRAS region (chr12:25,000,000-26,000,000):"
    echo ""

    for vcf in results/varscan/processed/*.snp.Somatic.hc.vcf; do
        sample=$(basename "$vcf" .snp.Somatic.hc.vcf)
        echo "--- $sample KRAS region ---"
        grep -v "^#" "$vcf" | awk '$1 == "12" || $1 == "chr12"' | \
            awk '$2 >= 25000000 && $2 <= 26000000 {print $1":"$2, $4">"$5}' | head -10
    done
fi

echo ""
echo "==================================================================="
echo "Checking broader gene regions for any somatic variants:"
echo "==================================================================="
echo ""

# KRAS: chr12:25,200,000-25,250,000
# MAP2K1: chr15:66,680,000-66,780,000
# TP53: chr17:7,571,000-7,590,000

echo "KRAS region (chr12:25.2M-25.25M) - any somatic SNPs:"
for vcf in results/varscan/processed/*.snp.Somatic.hc.vcf; do
    sample=$(basename "$vcf" .snp.Somatic.hc.vcf)
    count=$(grep -v "^#" "$vcf" | awk '($1 == "12" || $1 == "chr12") && $2 >= 25200000 && $2 <= 25250000' | wc -l)
    echo "  $sample: $count variants"
done

echo ""
echo "MAP2K1 region (chr15:66.68M-66.78M) - any somatic SNPs:"
for vcf in results/varscan/processed/*.snp.Somatic.hc.vcf; do
    sample=$(basename "$vcf" .snp.Somatic.hc.vcf)
    count=$(grep -v "^#" "$vcf" | awk '($1 == "15" || $1 == "chr15") && $2 >= 66680000 && $2 <= 66780000' | wc -l)
    echo "  $sample: $count variants"
done

echo ""
echo "TP53 region (chr17:7.571M-7.59M) - any somatic SNPs:"
for vcf in results/varscan/processed/*.snp.Somatic.hc.vcf; do
    sample=$(basename "$vcf" .snp.Somatic.hc.vcf)
    count=$(grep -v "^#" "$vcf" | awk '($1 == "17" || $1 == "chr17") && $2 >= 7571000 && $2 <= 7590000' | wc -l)
    echo "  $sample: $count variants"
done
