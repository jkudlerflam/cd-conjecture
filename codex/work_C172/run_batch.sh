#!/bin/bash
# C172 two-kink ED batch (local). Each line: tag, fam kwargs.
cd "$(dirname "$0")"
run(){ OMP_NUM_THREADS=1 python3 ed_kink.py "$2" 10 18 8 "$1" > "ed_$1.out" 2>&1; }
run P1 'a=0.9,Dz=2.0,g=0.3,hx=0.2,hz=0.05,Dx=0.03' &
run P2 'a=0.9,Dz=2.0,g=0.3,hx=0.2,hz=0.6,Dx=0.03' &
run P3 'a=0.9,Dz=2.0,g=0.3,hx=0.2,hz=1.0,Dx=0.03' &
run P4 'a=0.9,Dz=2.0,g=0.3,hx=0.2,hz=1.4,Dx=0.03' &
run P5 'a=0.9,Dz=2.0,g=0.6,G=0.3,hx=0.3,Dx=0.1,hz=0.05' &
run P6 'a=0.9,Dz=2.0,g=0.6,G=0.3,hx=0.3,Dx=0.1,hz=1.0' &
run P7 'a=0.9,Dz=2.0,hz=0.6' &
run P8 'a=0.9,Dz=2.0,g=0.35,gxy=-0.2,G=0.25,Gp=-0.15,Dx=0.2,Dy=0.1,hx=0.25,hy=-0.2,hz=0.05' &
wait
