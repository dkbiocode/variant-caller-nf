#!/usr/bin/env bash
sa() {
    sacct -X --format JobID,JobName,AllocCPUS,State,ExitCode,Elapsed,TimeLimit,Submit,Start,End $@
}
finished_jids=$(cut -f 3 reports/trace.txt  | grep -v native | paste -sd ',' -)
queued_jids=$(squeue -u $USER | grep amilan | awk '{print $1}' | paste -sd ',' - )
sa -j $finished_jids,$queued_jids
