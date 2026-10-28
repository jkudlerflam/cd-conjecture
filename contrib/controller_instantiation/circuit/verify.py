#!/usr/bin/env python3
"""Independent finite checks for the generated reversible gate stream."""
from __future__ import annotations
from collections import Counter, defaultdict
import itertools
import json
import math
from pathlib import Path
import random
import time

from sparse_circuit import (Graph, Localizer, branch_output, flag_word, logical_word,
                            make_circuit, projective_plane, reference_word, rotations,
                            simulate_packed)


def branch_labels(g):
    return [(0,a) for a in range(g.u)]+[(1,e) for e in range(g.n)]


def expected(g,occupation,branches):
    output=[];phase=0
    for right,(kind,label) in enumerate(branches):
        word,sign=branch_output(g,occupation,kind,label,bool(right))
        output+=word[:2*g.q]
        phase^=sign
        if kind:occupation^=1<<label
    return tuple(output)+logical_word(g,occupation),phase


def check_cases(g,cases,two_ports=True):
    circuit=make_circuit(g,two_ports)
    observed=simulate_packed(circuit,cases)
    for (occupation,branches),(word,sign,dirty) in zip(cases,observed):
        assert dirty==0, ('dirty scratch',occupation,branches)
        assert (word,sign)==expected(g,occupation,branches)
        assert Counter(word)==Counter(logical_word(g,occupation)+reference_word(g)*len(branches))
    return circuit,observed


def check_preparation(circuit):
    g=circuit.graph
    R=rotations(g)
    error=0.
    for stage in circuit.stages:
        state={0:1.}
        for op in stage.prep:
            name,*a=op
            new=defaultdict(complex)
            for word,amp in state.items():
                if name=='X':new[word^(1<<a[0])]+=amp
                elif name=='CX':new[word^((1<<a[1]) if word>>a[0]&1 else 0)]+=amp
                elif name=='CCX':new[word^((1<<a[2]) if word>>a[0]&1 and word>>a[1]&1 else 0)]+=amp
                else:
                    if name=='H':target=a[0];active=True;c=s=math.sqrt(.5);had=True
                    elif name=='CH':control,target=a;active=bool(word>>control&1);c=s=math.sqrt(.5);had=True
                    elif name=='R':index,control,target=a;active=bool(word>>control&1);c,s=(math.sqrt(float(x)) for x in R[index]);had=False
                    else:raise AssertionError(op)
                    if not active:new[word]+=amp;continue
                    if word>>target&1:
                        new[word]+=amp*(-c if had else c)
                        new[word^(1<<target)]+=amp*(s if had else -s)
                    else:
                        new[word]+=amp*c
                        new[word^(1<<target)]+=amp*s
            state={word:amp for word,amp in new.items() if abs(amp)>1e-14}
        support={}
        allowed_mask=(1<<stage.kind)|sum(1<<w for w in stage.label)
        for word,amp in state.items():
            assert word&~allowed_mask==0,'Branch-preparation scratch was not cleaned'
            kind=word>>stage.kind&1
            label=sum(((word>>w)&1)<<j for j,w in enumerate(stage.label))
            support[kind,label]=amp
        assert set(support)==set(branch_labels(g))
        for (kind,label),amp in support.items():
            error=max(error,abs(amp-1/math.sqrt(2*(g.n if kind else g.u))))
    return error


def check_compressions(g):
    labels=branch_labels(g)
    cases=[(occupation,(left,right)) for occupation in range(1<<g.n)
           for left in labels for right in labels]
    circuit,observed=check_cases(g,cases)
    supports=[{} for _ in range(1<<g.n)]
    lookup={}
    for (occupation,branches),(word,phase,dirty) in zip(cases,observed):
        amp=(-1 if phase else 1)
        for kind,label in branches:amp/=math.sqrt(2*(g.n if kind else g.u))
        assert word not in lookup,'Basis transducer was not injective'
        supports[occupation][word]=amp;lookup[word]=(occupation,amp)
    error=0.
    for occupation,support in enumerate(supports):
        error=max(error,abs(sum(abs(a)**2 for a in support.values())-1))
    epsilon=1/(2*math.sqrt(g.u*g.n))
    for side,site in [('left',0),('right',2*g.q)]:
        scale=epsilon*(1-1/(2*g.n) if side=='right' else 1)
        for edge,(upper,lower) in enumerate(g.edges):
            matrix=defaultdict(float)
            for occupation,support in enumerate(supports):
                for word,amp in support.items():
                    if word[site]!=g.u+lower:continue
                    changed=word[:site]+(upper,)+word[site+1:]
                    if changed in lookup:
                        row,left_amp=lookup[changed]
                        matrix[row,occupation]+=left_amp*amp
            for row in range(1<<g.n):
                for col in range(1<<g.n):
                    target=0.
                    if not(col>>edge&1) and row==col^(1<<edge):
                        phase=(col&((1<<edge)-1)).bit_count()
                        if side=='right':phase+=col.bit_count()
                        target=scale*(-1 if phase%2 else 1)
                    error=max(error,abs(matrix[row,col]-target))
    return circuit,{'basis_branch_cases':len(cases),'max_compression_error':error,
                    'max_preparation_error':check_preparation(circuit)}


def species(symbol):
    return symbol[1] if symbol[0]=='D' else symbol[3] if symbol[0]=='T' else None


def charge_pair(pair,q):
    return tuple(sum(species(x)==s for x in pair) for s in range(q))


def local_permutation(opcode,pair,g):
    a,b=pair
    if opcode=='SWAP':return b,a
    if opcode=='TYPE':
        if a[0]=='B':a=('H',0,a[1])
        elif a[0]=='H' and a[1]==0:a=('B',a[2])
    elif opcode=='INIT':
        if a==('blank',):a=('B',0)
        elif a==('B',0):a=('blank',)
    elif opcode in {'SHIFT','SHIFT_DAG'}:
        delta=1 if opcode=='SHIFT' else -1
        if a[0] in {'H','T'}:a=(a[0],(a[1]+delta)%g.m,*a[2:])
    elif opcode=='QUERY':
        if a[0]=='D' and b[0]=='H' and b[1]<g.species_width:
            b=('H',b[1],b[2]^((a[1]>>b[1])&1))
    elif opcode=='LOAD':
        if a[0]=='D' and b[0]=='H':a,b=('blank',),('T',b[1],b[2],a[1])
        elif a==('blank',) and b[0]=='T':a,b=('D',b[3]),('H',b[1],b[2])
    elif opcode=='EXCHANGE':
        if a[0]=='T' and a[2]==1 and b[0]=='D':a,b=('T',a[1],1,b[1]),('D',a[3])
    elif opcode=='X':
        if a[0]=='B':a=('B',a[1]^1)
    elif opcode=='CNOT':
        if a==('B',1) and b[0]=='B':b=('B',b[1]^1)
    else:raise ValueError(opcode)
    return a,b


def check_local_alphabet(g):
    sigma=[('D',s) for s in range(g.q)]+[('blank',),('B',0),('B',1)]
    sigma += [('H',j,c) for j in range(g.m) for c in (0,1)]
    sigma += [('T',j,c,s) for j in range(g.m) for c in (0,1) for s in range(g.q)]
    assert len(sigma)==g.sigma_dimension
    count=0
    for opcode in ('SWAP','TYPE','INIT','SHIFT','SHIFT_DAG','QUERY','LOAD','EXCHANGE','X','CNOT'):
        outputs=set()
        for pair in itertools.product(sigma,repeat=2):
            out=local_permutation(opcode,pair,g)
            assert charge_pair(pair,g.q)==charge_pair(out,g.q)
            outputs.add(out);count+=1
        assert len(outputs)==len(sigma)**2,('not a full-space permutation',opcode)
    return count


def check_localizer(circuit):
    """Actual nearest-neighbor permutation execution, including dirty idle records."""
    g=circuit.graph;local=Localizer(circuit);checks=0
    rng=random.Random(91317)
    # Use real data/neutral positions. Other registers carry arbitrary symbols;
    # route-and-restore must leave each of these exactly unchanged.
    q=circuit.data_sites
    for _ in range(20):
        control=q+rng.randrange(circuit.wires)
        left,right=rng.sample(range(q),2)
        for bit in (0,1):
            state=[('D',rng.randrange(g.q)) for _ in range(local.size)]
            state[control]=('B',bit)
            original=state.copy();wanted=state.copy()
            if bit:wanted[left],wanted[right]=wanted[right],wanted[left]
            for opcode,site in local.controlled_data_swap(control,left,right):
                state[site],state[site+1]=local_permutation(opcode,(state[site],state[site+1]),g)
            assert state==wanted
            checks+=1
    for data in (0,q-1):
        for target in (q,local.size-1):
            for s in range(g.q):
                for index in range(g.species_width):
                    state=[('D',rng.randrange(g.q)) for _ in range(local.size)]
                    state[data]=('D',s);state[target]=('B',1)
                    wanted=state.copy();wanted[target]=('B',1^((s>>index)&1))
                    for opcode,site in local.indexed_query(data,index,target):
                        state[site],state[site+1]=local_permutation(opcode,(state[site],state[site+1]),g)
                    assert state==wanted
                    checks+=1
    # Ordered pair gathering in both operand orders, including last site.
    for left,right in [(0,local.size-1),(local.size-1,0),(3,2),(2,3)]:
        for a,b in itertools.product((0,1),repeat=2):
            state=[('D',rng.randrange(g.q)) for _ in range(local.size)]
            state[left]=('B',a);state[right]=('B',b)
            wanted=state.copy();wanted[right]=('B',b^a)
            for opcode,site in local.pair('CNOT',left,right):
                state[site],state[site+1]=local_permutation(opcode,(state[site],state[site+1]),g)
            assert state==wanted
            checks+=1
    return checks


def check_toffoli():
    # Independent full coherent 3-bit unitary calculation for the five-gate
    # controlled-sqrt-X decomposition used by Localizer.abstract.
    V=(((1+1j)/2,(1-1j)/2),((1-1j)/2,(1+1j)/2))
    error=0.
    for initial in range(8):
        state={initial:1+0j}
        for kind,control,target in [('V',1,2),('X',0,1),('VD',1,2),('X',0,1),('V',0,2)]:
            new=defaultdict(complex)
            for word,amp in state.items():
                if not(word>>control&1):new[word]+=amp;continue
                if kind=='X':new[word^(1<<target)]+=amp;continue
                old=word>>target&1
                for bit in (0,1):
                    coefficient=V[bit][old]
                    if kind=='VD':coefficient=coefficient.conjugate()
                    new[(word&~(1<<target))|(bit<<target)]+=coefficient*amp
            state=new
        expected=initial^((1<<2) if initial&3==3 else 0)
        for word in range(8):error=max(error,abs(state.get(word,0)-(1 if word==expected else 0)))
    return error


def check_full_local_stream():
    """Execute the entire nearest-neighbor unitary, including coherent branches.

    This tests the actual localizer, neutral gates, initialized BLANK records,
    read/query machinery, routed charged swaps, and both inverse cleanups.
    """
    g=Graph(1,1,((0,0),));circuit=make_circuit(g);local=Localizer(circuit)
    sigma=[('D',s) for s in range(g.q)]+[('blank',),('B',0),('B',1)]
    sigma += [('H',j,c) for j in range(g.m) for c in (0,1)]
    sigma += [('T',j,c,s) for j in range(g.m) for c in (0,1) for s in range(g.q)]
    index={x:i for i,x in enumerate(sigma)}
    site_bits=(len(sigma)-1).bit_length();site_mask=(1<<site_bits)-1
    R=rotations(g)
    V=(((1+1j)/2,(1-1j)/2),((1-1j)/2,(1+1j)/2))
    def encode(values):return sum(index[x]<<(site_bits*i) for i,x in enumerate(values))
    def transitions(opcode,a,b):
        if opcode in {'H','Z'} and a[0]=='B':
            if opcode=='Z':return [(a,b,(-1 if a[1] else 1))]
            return [(('B',new),b,(-1 if a[1] and new else 1)/math.sqrt(2)) for new in (0,1)]
        if opcode in {'CV','CV_DAG','CH'} and a==('B',1) and b[0]=='B':
            if opcode=='CH':return [(a,('B',new),(-1 if b[1] and new else 1)/math.sqrt(2)) for new in (0,1)]
            return [(a,('B',new),V[new][b[1]].conjugate() if opcode=='CV_DAG' else V[new][b[1]]) for new in (0,1)]
        if opcode in {'ROT','ROT_DAG'} and a[0]=='H' and b==('B',1):
            cc,ss=(math.sqrt(float(x)) for x in R[a[1]])
            matrix=((cc,-ss),(ss,cc))
            if opcode=='ROT_DAG':matrix=((cc,ss),(-ss,cc))
            return [(('H',a[1],new),b,matrix[new][a[2]]) for new in (0,1)]
        if opcode in {'H','Z','CV','CV_DAG','CH','ROT','ROT_DAG'}:return [(a,b,1)]
        aa,bb=local_permutation(opcode,(a,b),g)
        return [(aa,bb,1)]
    # A finite pair transition cache prevents the state-vector check from
    # conflating gate semantics with the Boolean simulator.
    cache={}
    stream=list(local.stream())
    assert Counter(opcode for opcode,site in stream)==local.gate_counts()
    errors=[];max_support=0
    for occupation in range(2):
        original=[('D',s) for s in logical_word(g,occupation)+reference_word(g)*2]
        original += [('blank',)]*circuit.wires
        state={encode(original):1+0j}
        for opcode,site in stream:
            shift=site_bits*site;pair_mask=site_mask|(site_mask<<site_bits)
            new=defaultdict(complex)
            for word,amp in state.items():
                pair=(word>>shift)&pair_mask
                key=opcode,pair
                if key not in cache:
                    a=sigma[pair&site_mask];b=sigma[(pair>>site_bits)&site_mask]
                    cache[key]=[(index[aa]|(index[bb]<<site_bits),coefficient)
                                for aa,bb,coefficient in transitions(opcode,a,b) if coefficient]
                base=word&~(pair_mask<<shift)
                for output,coefficient in cache[key]:new[base|(output<<shift)]+=amp*coefficient
            state={word:amp for word,amp in new.items() if abs(amp)>1e-13}
            max_support=max(max_support,len(state))
        wanted={}
        for branches in itertools.product(branch_labels(g),repeat=2):
            word,phase=expected(g,occupation,branches)
            output=[('D',s) for s in word]+[('blank',)]*circuit.wires
            wanted[encode(output)]=(-1 if phase else 1)/2
        assert set(state)==set(wanted)
        errors.append(max(abs(state[word]-amp) for word,amp in wanted.items()))
    assert max(errors)<1e-12
    return {'local_gates':len(stream),'sites':local.size,
            'basis_columns_checked':2,'max_intermediate_support':max_support,
            'max_amplitude_error':max(errors),
            'all_scratch_returned_to_blank':True}


def main():
    start=time.monotonic()
    results={'status':'finite circuit checks only; not a counterexample'}
    graphs=[Graph(1,1,((0,0),)),Graph(2,2,((0,0),(0,1),(1,0))),
            Graph(2,2,((0,0),(0,1),(1,0),(1,1))),Graph(3,2,((0,0),(2,1)))]
    results['exhaustive_graphs']=[]
    for g in graphs:
        circuit,result=check_compressions(g)
        assert result['max_compression_error']<2e-14
        assert result['max_preparation_error']<2e-14
        results['exhaustive_graphs'].append(dict(u=g.u,v=g.v,N=g.n,**result))
    tiny=make_circuit(graphs[0])
    results['local_full_space_permutation_inputs']=check_local_alphabet(graphs[0])
    results['local_routing_checks']=check_localizer(tiny)
    results['toffoli_max_error']=check_toffoli();assert results['toffoli_max_error']<1e-14
    results['full_nearest_neighbor_unitary']=check_full_local_stream()
    g=projective_plane(2)
    assert g.u==g.v==7 and g.n==21
    assert set(Counter(a for a,b in g.edges).values())=={3}
    assert set(Counter(b for a,b in g.edges).values())=={3}
    assert all(len({v for a,v in g.edges if a==x}&{v for a,v in g.edges if a==y})==1
               for x in range(g.u) for y in range(x))
    rng=random.Random(20261028);labels=branch_labels(g)
    cases=[(rng.randrange(1<<g.n),(rng.choice(labels),rng.choice(labels))) for _ in range(1024)]
    cases += [(occupation,(left,right)) for occupation in (0,(1<<g.n)-1)
              for left in labels for right in labels]
    circuit,observed=check_cases(g,cases)
    results['projective_plane_p2']={'basis_branch_cases':len(cases),
                                  'max_preparation_error':check_preparation(circuit),
                                  'all_species_counts_preserved':True,
                                  'all_neutral_scratch_cleared':True}
    assert results['projective_plane_p2']['max_preparation_error']<1e-14
    results['elapsed_seconds']=round(time.monotonic()-start,3)
    here=Path(__file__).resolve().parent
    (here/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    (here/'p2_manifest.json').write_text(json.dumps(circuit.manifest(),indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
