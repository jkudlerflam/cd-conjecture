#!/bin/bash
# production set, 2026-11-18 (resumable; re-running this only restarts missing stages)
cd ~/cc_spinchain/v3
L="chis=32,48,64,96,128,160,200,256,320 maxh=6"
for h in 1.738 1.740 1.741 1.742 1.743 1.744 1.745 1.746 1.747 1.748 1.749 1.750 1.752 1.755 1.758; do sbatch -J P_$h run3.sbatch model=target h=$h $L branch=pin >/dev/null; done
sub(){ n=$1; shift; sbatch -J $n run3.sbatch model=ctrl $L branch=pin "$@" >/dev/null; }
sub C1 tA=1 hA=0 J=1 g=1; sub C2 tA=1 hA=0 J=0.35 g=0.35; sub C3 tA=1 hA=0 J=2.5 g=2.5; sub C4 tA=1 JzA=-0.39 hA=0 J=1 g=1
sub C5 tA=1 hA=0 J=1 g=0.8; sub C6 tA=1 hA=0 J=1 g=0.97; sub C7 tA=0 hA=0 J=1 g=1; sub C8 tA=1 hA=0 J=0 g=1
for h in 1.740 1.744 1.746 1.748 1.750; do sbatch -J U_$h run3.sbatch model=target h=$h $L branch=cplx >/dev/null; done
for g in 0.98 0.99 1.01 1.02; do sub Cg$g tA=1 hA=0 J=1 g=$g; done
DN=$(python3 -c "print(\",\".join(\"%.4f\"%(1.754-0.001*i) for i in range(17)))"); UP=$(python3 -c "print(\",\".join(\"%.4f\"%(1.738+0.001*i) for i in range(17)))")
for chi in 64 128; do
sbatch -J HD_$chi run3.sbatch model=target hs=$DN chis=$chi chi0=32 branch=pin tag=_hdn$chi maxh=1.5 corr=0 nmix=0 >/dev/null
sbatch -J HU_$chi run3.sbatch model=target hs=$UP chis=$chi chi0=32 branch=cplx seed=1e-4 tag=_hup$chi maxh=1.5 corr=0 nmix=0 >/dev/null
done
for h in 1.739 1.7395 1.740 1.7405 1.741 1.742 1.744; do sbatch -J B_$h run3.sbatch model=target h=$h $L branch=pinB >/dev/null; done
for h in 1.7455 1.7465 1.7475; do sbatch -J P_$h run3.sbatch model=target h=$h $L branch=pin >/dev/null; done
