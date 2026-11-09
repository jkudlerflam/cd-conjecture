export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
for pv in "0.5 -2.4" "0.5 -2.2" "0.5 -2.0" "0.5 -1.8" "0.35 -2.0" "0.35 -1.8"; do
  set -- $pv
  p=$(python3 -c "from model import spin_from_fermion as s; print(','.join('%.10f'%x for x in s(1.0,1.0,$1,$2,$2)))")
  for L in 32 64; do python3 dmrg_pt.py $L 16,32,64,96,128,128 local_out/dmrg_$L.jsonl ph1.00_D$1_V$2 spin:$p; done
done
