#!/bin/bash
# usage: body.sh <task list> <parallel tasks>
cd /scratch/lustre/jkudlerflam/c187d || exit 1
echo "start $1 on $(hostname) at $(date)"
cat "$1" | xargs -P "$2" -I CMD sh -c 'CMD'
echo "payload finished $1 at $(date)"
