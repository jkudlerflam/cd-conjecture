#!/bin/bash
# Open-chain finite DMRG, second method on the c ~ 2 region (2026-11-19).  Run from ~/cc_spinchain/v3.
CH=chis=64,128,256,384
CHL=chis=64,128,256,384,512
sub() { sbatch -J "$1" "${@:2}"; }
for pt in "T1 model=target d=0.1875 h=1.783" "T2 model=target d=0.125 h=1.7405" "T3 model=target d=0.1875 h=1.772"; do
  set -- $pt; n=$1; shift
  sub o${n}_s runo.sbatch "$@" L=32,64,96,128,192,256 $CH
  sub o${n}_384 runo.sbatch "$@" L=384 $CHL
  sub o${n}_512 runo.sbatch "$@" L=512 $CHL
done
for pt in "C2 model=ctrl tA=1 hA=0.3 tB=1 J=0 g=0.3" "C15 model=ctrl tA=1 hA=0.3 J=1 g=1" "C1 model=ctrl tA=1 hA=0.3 J=1 g=3"; do
  set -- $pt; n=$1; shift
  sub o${n}_s runo.sbatch "$@" L=32,64,96,128,192 $CHL
  sub o${n}_256 runo.sbatch "$@" L=256 $CHL
  sub o${n}_384 runo.sbatch "$@" L=384 $CHL
done
# --- chi extension (2026-11-20): S(L/2) not converged at chi<=512 in the c~2 region and for XXxXX
X=chis=64,128,256,384,512,768,1024
for pt in "T1 model=target d=0.1875 h=1.783" "T2 model=target d=0.125 h=1.7405"; do
  set -- $pt; n=$1; shift
  sub o${n}_x1 runob.sbatch "$@" L=96,128 $X
  sub o${n}_x2 runob.sbatch "$@" L=192 $X
  sub o${n}_x3 runob.sbatch "$@" L=256 $X
done
sub oC2_x1 runob.sbatch model=ctrl tA=1 hA=0.3 tB=1 J=0 g=0.3 L=96,128 $X
sub oC2_x2 runob.sbatch model=ctrl tA=1 hA=0.3 tB=1 J=0 g=0.3 L=192 $X
