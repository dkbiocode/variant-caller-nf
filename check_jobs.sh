#!/usr/bin/env bash

# Configurable refresh interval (default 2 seconds, or pass as argument)
REFRESH_INTERVAL=${1:-2}

# Function to get sacct output
sa() {
    sacct -X --format JobID,JobName%30,AllocCPUS,State,ExitCode,Elapsed,TimeLimit,Submit,Start,End "$@"
}

# Trap Ctrl-C to exit cleanly
trap 'tput cnorm; exit 0' INT TERM

# Hide cursor for cleaner display
tput civis

while true; do
    # Move cursor to top-left without clearing (no flicker)
    tput cup 0 0

    # Header with timestamp
    echo "======================================================================"
    echo "  SLURM Job Monitor - $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  Refresh: ${REFRESH_INTERVAL}s"
    echo "======================================================================"
    echo ""

    # Get finished job IDs from Nextflow trace
    if [[ -f reports/trace.txt ]]; then
        finished_jids=$(cut -f 3 reports/trace.txt 2>/dev/null | grep -v native | grep -v '^$' | paste -sd ',' -)
    else
        finished_jids=""
    fi

    # Get queued job IDs from squeue
    queued_jids=$(squeue -u "$USER" 2>/dev/null | grep amilan | awk '{print $1}' | paste -sd ',' -)

    # Combine job IDs
    if [[ -n "$finished_jids" && -n "$queued_jids" ]]; then
        all_jids="$finished_jids,$queued_jids"
    elif [[ -n "$finished_jids" ]]; then
        all_jids="$finished_jids"
    elif [[ -n "$queued_jids" ]]; then
        all_jids="$queued_jids"
    else
        echo "No jobs found"
        echo ""
        echo "Press Ctrl-C to exit"
        sleep "$REFRESH_INTERVAL"
        continue
    fi

    # ANSI color codes
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[0;33m'
    CYAN='\033[0;36m'
    BLUE='\033[0;34m'
    MAGENTA='\033[0;35m'
    NC='\033[0m' # No Color

    # Show job status with color coding
    echo "Jobs Status:"
    echo "----------------------------------------------------------------------"

    # Get sacct output and colorize only the State column
    sa -j "$all_jids" 2>/dev/null | while IFS= read -r line; do
        # Skip empty lines
        [[ -z "$line" ]] && continue

        # Color code only the state word, not the entire line
        # Use printf instead of echo to preserve formatting and enable colors
        printf '%s\n' "$line" | sed -E \
            -e "s/RUNNING/$(printf "${GREEN}")RUNNING$(printf "${NC}")/g" \
            -e "s/PENDING/$(printf "${YELLOW}")PENDING$(printf "${NC}")/g" \
            -e "s/COMPLETED/$(printf "${CYAN}")COMPLETED$(printf "${NC}")/g" \
            -e "s/FAILED/$(printf "${RED}")FAILED$(printf "${NC}")/g" \
            -e "s/TIMEOUT/$(printf "${MAGENTA}")TIMEOUT$(printf "${NC}")/g" \
            -e "s/CANCELLED/$(printf "${BLUE}")CANCELLED$(printf "${NC}")/g"
    done

    echo ""

    # Summary statistics with color
    echo "----------------------------------------------------------------------"
    echo "Summary:"
    if [[ -n "$queued_jids" ]]; then
        queued_count=$(echo "$queued_jids" | tr ',' '\n' | wc -l)
        echo -e "  ${YELLOW}Queued/Running: $queued_count${NC}"
    else
        echo "  Queued/Running: 0"
    fi

    if [[ -n "$finished_jids" ]]; then
        finished_count=$(echo "$finished_jids" | tr ',' '\n' | wc -l)

        # Count by status - ensure we get a clean integer
        completed=$(sa -j "$all_jids" 2>/dev/null | grep -c "COMPLETED" | tr -d '\n' || echo 0)
        failed=$(sa -j "$all_jids" 2>/dev/null | grep -c "FAILED" | tr -d '\n' || echo 0)

        # Remove any whitespace and ensure it's a number
        completed=${completed// /}
        failed=${failed// /}
        : ${completed:=0}
        : ${failed:=0}

        echo -e "  ${CYAN}Completed: $completed${NC}"
        if [[ "$failed" -gt 0 ]] 2>/dev/null; then
            echo -e "  ${RED}Failed: $failed${NC}"
        fi
        echo "  Total Finished: $finished_count"
    else
        echo "  Finished: 0"
    fi

    echo ""
    echo "Press Ctrl-C to exit | Updating every ${REFRESH_INTERVAL} seconds..."
    echo ""
    echo -e "Legend: ${GREEN}RUNNING${NC} | ${YELLOW}PENDING${NC} | ${CYAN}COMPLETED${NC} | ${RED}FAILED${NC} | ${MAGENTA}TIMEOUT${NC} | ${BLUE}CANCELLED${NC}"

    # Pad with blank lines to clear any leftover text from previous iteration
    for i in {1..20}; do echo ""; done

    # Wait before next update
    sleep "$REFRESH_INTERVAL"
done
