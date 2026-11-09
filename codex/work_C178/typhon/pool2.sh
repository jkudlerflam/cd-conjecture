#!/bin/bash
# second worker pool inside job 17603981 (srun --overlap), taking medium.txt from the end; locks prevent duplicates
cd /scratch/lustre/jkudlerflam/c178
tac medium.txt | sed 's/^/4 /' | xargs -P 6 -L 1 ./run_task.sh
echo pool2 done $(date) >> logs/done.txt
