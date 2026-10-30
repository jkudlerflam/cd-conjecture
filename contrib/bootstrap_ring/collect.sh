#!/bin/bash
# Assemble the final result set from the campaign directories (run on Typhon in the work dir).
# out/main : first campaign (ED, local, s=1, ring wraparound); its s=2 runs used SCS and are superseded
# v3/out   : d=3 runs with the energy-density operator set
# v4/out   : regularized-Clarabel s=2 runs and ingredient tests (supersede earlier s=2 runs)
set -e
rm -rf final && mkdir -p final
for f in out/main/*.json; do case "$f" in *_s2.json) ;; *) ln -s "../$f" final/ ;; esac; done
for f in v3/out/SU3_swap*_k5_* v3/out/SU3_swap*_s2_oh.json v3/out/SU3_swap*_honly.json v3/out/random_real_d3*_oh.json; do ln -sf "../$f" final/; done
for f in v4/out/*.json; do ln -sf "../$f" final/; done
