#!/bin/bash
# Watchdog for the full-tissue CytoSignal run on a RAM/disk-constrained box.
# Kills the run if free disk < MIN_DISK_GB or swap used > MAX_SWAP_GB, logging the trajectory.
# Usage: watchdog.sh <log_file>
LOG="${1:-/tmp/cyto_watchdog.log}"
MIN_DISK_GB=3
MAX_SWAP_GB=20
PAT='run_cytosignal.R'
echo "$(date +%H:%M:%S) watchdog start (min_disk=${MIN_DISK_GB}G max_swap=${MAX_SWAP_GB}G)" > "$LOG"
for i in $(seq 1 120); do pgrep -f "$PAT" >/dev/null && break; sleep 1; done   # wait for run to appear
while true; do
  pid=$(pgrep -f "$PAT" | head -1)
  if [ -z "$pid" ]; then echo "$(date +%H:%M:%S) run process gone; watchdog exit" >> "$LOG"; break; fi
  free_gb=$(( $(df -k /Users/jiayifan | tail -1 | awk '{print $4}') / 1024 / 1024 ))
  swap_mb=$(sysctl -n vm.swapusage | sed -E 's/.*used = ([0-9.]+)([MG]).*/\1 \2/' | awk '{printf "%.0f", ($2=="G")?$1*1024:$1}')
  swap_gb=$(( swap_mb / 1024 ))
  rss_gb=$(ps -o rss= -p "$pid" 2>/dev/null | awk '{printf "%.1f", $1/1024/1024}')
  pressure=$(memory_pressure 2>/dev/null | awk -F': ' '/System-wide memory free percentage/{print $2}')
  echo "$(date +%H:%M:%S) pid=$pid rss=${rss_gb}G free_disk=${free_gb}G swap_used=${swap_gb}G memfree=${pressure}" >> "$LOG"
  if [ "$free_gb" -lt "$MIN_DISK_GB" ]; then echo "$(date +%H:%M:%S) *** BAIL: free disk ${free_gb}G < ${MIN_DISK_GB}G — killing $pid" >> "$LOG"; kill -9 "$pid" 2>/dev/null; pkill -9 -f "$PAT"; break; fi
  if [ "$swap_gb" -gt "$MAX_SWAP_GB" ]; then echo "$(date +%H:%M:%S) *** BAIL: swap ${swap_gb}G > ${MAX_SWAP_GB}G — killing $pid" >> "$LOG"; kill -9 "$pid" 2>/dev/null; pkill -9 -f "$PAT"; break; fi
  sleep 15
done
echo "$(date +%H:%M:%S) watchdog done" >> "$LOG"
