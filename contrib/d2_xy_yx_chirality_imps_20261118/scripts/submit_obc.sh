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
# --- 2026-11-23: odd-even energies across the lower edge of the wedge
for h in 1.734 1.737 1.7405 1.744; do for p in 0 1; do sub oP${p}_8_$h runo.sbatch model=target d=0.125 h=$h par=$p L=64,96,128,192 chis=64,128,256,512; done; done
for h in 1.774 1.778 1.783; do for p in 0 1; do sub oP${p}_3_$h runo.sbatch model=target d=0.1875 h=$h par=$p L=64,96,128,192 chis=64,128,256,512; done; done
# 2026-11-23b: open-chain replacement for the failed D1_/D2_ iMPS scans (iDMRG fell into the reversed state)
for h in 1.680 1.685 1.690 1.695 1.700 1.705 1.710 1.715 1.720 1.725 1.730; do sub oD1_$h runo.sbatch model=target d=0.0625 h=$h par=0 L=64,128 chis=64,128,256; done
for h in 1.705 1.710 1.715 1.720 1.725 1.730 1.735 1.740; do sub oD2_$h runo.sbatch model=target d=0.09375 h=$h par=0 L=64,128 chis=64,128,256; done
# 2026-11-23c: fine scans at the onset of the gapless region below delta = 1/8
for h in 1.727 1.729 1.731 1.733 1.735; do sub oF2_$h runo.sbatch model=target d=0.09375 h=$h par=0 L=64,96,128,192 chis=64,128,256,512; done
for h in 1.733 1.736 1.739 1.742 1.745; do sub oF1_$h runo.sbatch model=target d=0.0625 h=$h par=0 L=64,96,128,192 chis=64,128,256,512; done
