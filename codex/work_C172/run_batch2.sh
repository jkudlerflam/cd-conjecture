#!/bin/bash
# C172 second ED batch: order-one charge-2 amplitude at Dz = 2, and Dz = 1.6 with its parent control.
cd "$(dirname "$0")"
run(){ OMP_NUM_THREADS=1 python3 ed_kink.py "$2" 10 18 8 "$1" > "ed_$1.out" 2>&1; }
run Q1 'a=0.9,Dz=2.0,g=0.9,G=0.4,hx=0.3,Dx=0.1,hz=0.05' &
run Q2 'a=0.9,Dz=2.0,g=0.8,gxy=0.4,G=0.5,Gp=-0.3,Dy=0.3,hx=0.4,hy=0.2,hz=0.3' &
run R0 'a=0.9,Dz=1.6' &
run R1 'a=0.9,Dz=1.6,g=0.6,G=0.3,hx=0.3,Dx=0.1,hz=0.05' &
run R2 'a=0.9,Dz=1.6,g=0.35,gxy=-0.2,G=0.25,Gp=-0.15,Dx=0.2,Dy=0.1,hx=0.25,hy=-0.2,hz=0.05' &
run S0 'a=2.2,Dz=2.0' &
run S1 'a=2.2,Dz=2.0,g=0.7,G=0.3,hx=0.3,Dx=0.1,hz=0.05' &
wait
