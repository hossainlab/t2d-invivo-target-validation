#!/usr/bin/env bash
# run_integration_chain.sh — steps that combine bulk candidates with scRNA, run after the bulk chain
# (01-05, 12a, 11) and the scRNA chain (06a-09) have finished.
#   10 integration (liver, kidney) -> 13 co-localization / correlation (liver, kidney) -> 12 prioritization
# Usage (from project root): bash scripts/run_integration_chain.sh
# Logs: logs/<step>_<tissue>.log. Stops at the first failing step.
# Do not edit any of these scripts while the chain is running (Rscript reads scripts incrementally).

set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd "$(dirname "$0")/.." || exit 1

run() {  # run <log_name> <script> [args...]
  local log="logs/$1.log"; shift
  echo "[$(date +%H:%M:%S)] START $*"
  "$RS" "$@" > "$log" 2>&1
  local code=$?
  echo "[$(date +%H:%M:%S)] EXIT $code  $*  (log: $log)"
  if [ $code -ne 0 ]; then tail -25 "$log"; exit $code; fi
}

for t in liver kidney; do run "10_integration_$t"    scripts/10_integration_bulk_sc.R  "$t"; done
for t in liver kidney; do run "13_colocalization_$t" scripts/13_sc_colocalization.R   "$t"; done
run 12_candidate_prioritization scripts/12_candidate_prioritization.R
echo "[$(date +%H:%M:%S)] integration chain finished"
