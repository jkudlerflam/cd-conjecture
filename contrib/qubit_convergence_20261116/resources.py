#!/usr/bin/env python3
"""Explicit resource model from real-gauge arrays and archived sweep times."""
import json
from continuation import ROOT,HERE,storage,digest


def main():
    p=ROOT/'codex/work_C195/typhon_out_x2048/sym40_n316_L512_P1_rand1.json'
    r=json.loads(p.read_text()); times=[h['t'] for h in r['hist']]
    rows=[]
    for chi in (2048,2560,3072):
        ratio=chi/2048
        # Timing exponents 2 to 3 are an estimate, not a new large-chi benchmark.
        lo=min(times)*ratio**2/3600; hi=max(times)*ratio**3/3600
        rows.append(dict(chi=chi,**storage(512,chi),sweep_hours_range=[lo,hi],
                         stage_3_to_8_sweeps_hours=[3*lo,8*hi],
                         final_measurement_hours_cubic_estimate=r['t_meas']*ratio**3/3600,
                         uncertainty='allow a further factor two for matvec count, IO, hardware and measurements'))
    out=dict(host_observed=dict(memory_bytes=51539607552,memory_GiB=48,cpu_count=18,
        cuda_available=False,production_checkpoints_found=False,
        cluster_read_only_probe='ssh -o BatchMode=yes -o ConnectTimeout=5 typhon true: connection timed out',
        date='2026-11-16 America/New_York'),source=str(p.relative_to(ROOT)),source_sha256=digest(p),
        observed_chi2048_sweep_seconds=times,rows=rows,
        selected=dict(chi=2560,reason='next stage; lower memory than 3072',
            host_allocation_GB=256,gpu_recommendation='one CUDA GPU with 80 GB, benchmark before production',
            concurrent_runs=1,checkpoint_scratch_GB=300,
            both_parities_elapsed_hours_nominal=[34,66],
            both_parities_with_factor_two_margin_hours=[68,132],
            correlation_warning='C193 converts the complete MPS to complex128 on device: about 103 GB at chi2560 before stacks. An 80 GB GPU cannot run that unchanged; CPU or a separately validated streaming implementation is needed.'),
        production_status='NOT RUN: current host below memory floor, no CUDA or accessible remote checkpoints')
    (HERE/'resources.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))


if __name__=='__main__':main()
