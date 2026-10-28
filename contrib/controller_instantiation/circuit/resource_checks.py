#!/usr/bin/env python3
"""Check conservative source-loop resource bounds, without huge local streams."""
import gc
import json
from pathlib import Path
import time

from sparse_circuit import Graph, Localizer, make_circuit, projective_plane


def bounds(u, v, edges):
    if min(u, v, edges) < 1 or edges > u*v:
        raise ValueError('Require a nonempty finite simple bipartite graph')
    q = u+v
    L = edges+4*q
    w = max(1, (q-1).bit_length(), (max(u,edges)-1).bit_length())
    A = L*L*w
    # These are the unsimplified totals in RESOURCE_BOUND.md, with one block
    # for each of the two ports. They are not fitted to the generated counts.
    wire_loop_bound = 2*(9*A+22*L*w+2*L*w*w+w+1)
    gate_loop_bound = 84*A+L*L+12*L*w*w+26*L*w+2*w+6
    assert wire_loop_bound <= 32*A
    assert gate_loop_bound <= 128*A
    return dict(L=L,w=w,A=A,neutral_bits=32*A,abstract_gates=128*A,
                total_data_sites=33*A,local_gates=70000*A*A,
                cursor_program_symbols=2310000*A**3,
                preloaded_tape_sites=2310035*A**3,
                dispatcher_cycle_steps=30000000000000*A**6,
                unsimplified_wire_bound=wire_loop_bound,
                unsimplified_abstract_gate_bound=gate_loop_bound)


def main():
    start=time.perf_counter()
    cases=[('one-edge',Graph(1,1,((0,0),))),
           ('isolated-upper',Graph(8,1,((7,0),))),
           ('complete-3x2',Graph(3,2,tuple((a,b) for a in range(3) for b in range(2)))),
           ('projective-p2',projective_plane(2)),
           ('projective-p3',projective_plane(3))]
    results=[]
    for label,g in cases:
        c=make_circuit(g)
        manifest=c.manifest()
        upper=bounds(g.u,g.v,g.n)
        exact=dict(neutral_bits=c.wires,total_data_sites=c.data_sites+c.wires,
                   abstract_gates=manifest['abstract_gate_count'],
                   local_gates=manifest['local_gate_count'])
        for key,value in exact.items():assert value<=upper[key],(label,key,value,upper[key])
        assert c.wires<=upper['unsimplified_wire_bound']
        assert manifest['abstract_gate_count']<=upper['unsimplified_abstract_gate_bound']
        if label=='one-edge':
            cursor=0;program_size=0;local_size=0
            for opcode,left in Localizer(c).stream():
                program_size+=abs(left-cursor)+1
                cursor=left;local_size+=1
            assert local_size==exact['local_gates']
            assert program_size<=local_size*(exact['total_data_sites']-1)
            assert program_size<=upper['cursor_program_symbols']
            exact['cursor_program_symbols']=program_size
            exact['preloaded_tape_sites']=program_size+2+exact['total_data_sites']
        results.append(dict(case=label,q=g.q,N=g.n,exact=exact,bounds=upper))
        del c,manifest
        gc.collect()
    # Parameter-only checks do not construct large projective graphs/circuits.
    parameter_checks=0
    for p in range(2,1001):
        u=p*p+p+1;q=2*u;N=(p+1)*u
        upper=bounds(u,u,N)
        assert upper['L']<10*p**3
        assert upper['w']<=3*p
        assert upper['A']<=300*p**7
        assert upper['preloaded_tape_sites'] < 63000000000000*p**21
        assert upper['preloaded_tape_sites'] < (1<<50)*p**21
        parameter_checks+=1
    out=dict(status='Source-loop bounds for circuit plus preloaded dispatcher only',
             excluded='Autonomous program generation, verifier, history clock, endpoint transport, full core and wrapper',
             cases=results,projective_parameter_only_checks=parameter_checks,
             all_passed=True,runtime_seconds=time.perf_counter()-start)
    path=Path(__file__).resolve().with_name('resource_results.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
