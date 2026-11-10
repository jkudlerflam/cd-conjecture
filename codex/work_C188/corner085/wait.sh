#!/bin/bash
# Wait for both local run chains, then print the last lines of each log.
cd "$(dirname "$0")"
A=$(cat pid_A.txt); B=$(cat pid_B.txt)
while kill -0 "$A" 2>/dev/null || kill -0 "$B" 2>/dev/null; do sleep 60; done
for f in log_L48.txt log_L96.txt log_L192.txt; do echo "== $f"; tail -n 3 "$f" | cut -c1-200; done
