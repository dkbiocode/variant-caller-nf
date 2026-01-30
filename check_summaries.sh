#!/usr/bin/env bash
# Check the variant summary files

echo "==================================================================="
echo "Checking Variant Summary Files"
echo "==================================================================="
echo ""

if [[ -d "results/varscan/summaries" ]]; then
    echo "Summary files:"
    ls -lh results/varscan/summaries/
    echo ""

    echo "==================================================================="
    echo "Preview of each summary (first 10 lines):"
    echo "==================================================================="
    echo ""

    for tsv in results/varscan/summaries/*_variants.tsv; do
        if [[ -f "$tsv" ]]; then
            echo "--- $(basename "$tsv") ---"
            head -10 "$tsv"
            echo ""
        fi
    done
else
    echo "Summaries directory not found. The SUMMARIZE_VARIANTS process may not have run."
fi

echo "==================================================================="
echo "Checking for aggregated results:"
echo "==================================================================="
echo ""

if [[ -f "results/varscan/all_variants_summary.tsv" ]]; then
    echo "Found all_variants_summary.tsv"
    echo ""
    echo "Total lines:"
    wc -l results/varscan/all_variants_summary.tsv
    echo ""
    echo "Preview:"
    head -20 results/varscan/all_variants_summary.tsv
    echo ""
    echo "Summary stats:"
    echo "Total variants: $(tail -n +2 results/varscan/all_variants_summary.tsv | wc -l)"
    echo "Somatic variants: $(tail -n +2 results/varscan/all_variants_summary.tsv | grep Somatic | wc -l)"
else
    echo "all_variants_summary.tsv not found. The AGGREGATE_RESULTS process may not have run."
fi
