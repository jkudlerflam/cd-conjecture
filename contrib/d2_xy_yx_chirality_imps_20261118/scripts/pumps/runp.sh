#!/bin/bash
# D ladder for pumps.py: runp.sh <tag> <Dlist comma> <pumps args...>; each D warm-starts from the previous A.
tag=$1; Dl=$2; shift 2
mkdir -p out; prev=""
for D in ${Dl//,/ }; do
  if [ -z "$prev" ]; then sched="Dsched=$(python3 -c "print(','.join(str(x) for x in range(8,$D,8))+',$D' if $D>8 else '$D')")"; ld=""; else sched="Dsched=$D"; ld="load=$prev"; fi
  echo "=== $tag D=$D $(date)"
  python3 pumps.py "$@" $sched $ld save=out/${tag}_D$D > out/${tag}_D$D.txt 2>&1 || { echo "FAIL D=$D"; tail -5 out/${tag}_D$D.txt; exit 1; }
  grep "^D=$D " out/${tag}_D$D.txt
  prev=out/${tag}_D${D}_A.npy
done
echo "done $tag $(date)"
