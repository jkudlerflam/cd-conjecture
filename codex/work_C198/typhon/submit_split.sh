#!/bin/bash
# Submit the nine C198 runs as separate jobs (threads, memory per run).
cd /scratch/lustre/jkudlerflam/c19x/c198
while read tag t m; do sbatch --job-name=c198_$tag --cpus-per-task=$t --mem=$m typhon/c198_one.slurm $tag; done <<'LIST'
sym316_L1024 8 120G
sym316_L512 8 80G
c2_L1024 4 80G
sym316_L384 4 48G
c2_L512 4 48G
sym316_L256 4 32G
c1_L512 4 32G
sym316_L192 4 24G
sym316_L128 4 16G
LIST
