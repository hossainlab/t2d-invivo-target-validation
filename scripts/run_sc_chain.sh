#!/usr/bin/env bash
# run_sc_chain.sh — run the scRNA steps after 06a_sc_demux_embed.R has finished.
#   06b tissue assignment -> per tissue (liver, kidney): 07 QC/cluster -> 08 annotate -> 09 pseudobulk DE -> 10 bulk integration
# Usage (from project root): bash scripts/run_sc_chain.sh [start_step]
#   start_step: 06b (default), 07, 08, 09 or 10 — resume from that step.
# Logs: logs/<script>[_<tissue>].log. Stops at the first failing step.

set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
START="${1:-06b}"
cd "$(dirname "$0")/.." || exit 1

order=(06b 07 08 09 10)
started=0
should_run() { [ "$1" = "$START" ] && started=1; [ $started -eq 1 ]; }

run() {  # run <log_name> <script> [args...]
  local log="logs/$1.log"; shift
  echo "[$(date +%H:%M:%S)] START $*"
  "$RS" "$@" > "$log" 2>&1
  local code=$?
  echo "[$(date +%H:%M:%S)] EXIT $code  $*  (log: $log)"
  if [ $code -ne 0 ]; then tail -25 "$log"; exit $code; fi
}

for step in "${order[@]}"; do
  should_run "$step" || continue
  case "$step" in
    06b) run 06b_sc_tissue_assign scripts/06b_sc_tissue_assign.R ;;
    07)  for t in liver kidney; do run "07_sc_qc_cluster_$t"    scripts/07_sc_qc_cluster.R    "$t"; done ;;
    08)  for t in liver kidney; do run "08_sc_annotate_$t"      scripts/08_sc_annotate.R      "$t"; done ;;
    09)  for t in liver kidney; do run "09_sc_pseudobulk_DE_$t" scripts/09_sc_pseudobulk_DE.R "$t"; done ;;
    10)  for t in liver kidney; do run "10_integration_$t"      scripts/10_integration_bulk_sc.R "$t"; done ;;
  esac
done
echo "[$(date +%H:%M:%S)] scRNA chain finished"
