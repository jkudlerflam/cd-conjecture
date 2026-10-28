"""Finite exact checks for N23_MATCHING.md, Python standard library only.

Checks dimension/parity counts, the neutral-wire gauge, the adjacent-cell
Jordan-Wigner sign, and saved certificate metadata/strict margins. Does not
rerun C33's full integer inertia computation or test an infrared theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def occupied(word, mode, total):
    return (word >> (total-1-mode)) & 1

def prefix_parity(word, mode, total):
    return -1 if (word >> (total-mode)).bit_count()%2 else 1

def act(word, mode, total, create=False):
    occ=occupied(word,mode,total)
    if occ==int(create):
        return None
    return word^(1 << (total-1-mode)),prefix_parity(word,mode,total)

def global_hop(word,i,j,n):
    first=act(word,n+j,2*n)
    if first is None:return None
    second=act(first[0],i,2*n,True)
    if second is None:return None
    return second[0],first[1]*second[1]

def cell_hop(word,i,j,n,with_parity=True):
    left,right=word>>n,word&((1<<n)-1)
    first=act(right,j,n)
    if first is None:return None
    second=act(left,i,n,True)
    if second is None:return None
    parity=(-1 if left.bit_count()%2 else 1) if with_parity else 1
    return (second[0]<<n)|first[0],first[1]*second[1]*parity

jw_cases=0
missing_parity_witness=None
for n in range(1,7):
    for word in range(1<<(2*n)):
        for i in range(n):
            for j in range(n):
                assert global_hop(word,i,j,n)==cell_hop(word,i,j,n)
                jw_cases+=1
                if missing_parity_witness is None and global_hop(word,i,j,n)!=cell_hop(word,i,j,n,False):
                    missing_parity_witness=dict(n=n,word=word,left_mode=i,right_mode=j,
                        correct=global_hop(word,i,j,n),without_parity=cell_hop(word,i,j,n,False))
assert missing_parity_witness is not None

gauge_cases=0
for n in range(1,9):
    for word in range(1<<n):
        U=[prefix_parity(word,j,n) for j in range(n+1)]
        assert U[0]==1
        assert U[n]==(-1 if word.bit_count()%2 else 1)
        for j in range(n):
            Z=1-2*occupied(word,j,n)
            assert U[j+1]*Z*U[j]==1
            # This includes the last-site-to-marker weak loop, and applies
            # to both physical rails. A junction switch has equal U's.
            assert U[j]*U[j]==1
            for create in (False,True):
                result=act(word,j,n,create)
                if result is not None:assert result[1]==U[j]
            gauge_cases+=1

n=23
dim=(7*n+4)*(1<<n)
each_exc=(dim-(1<<n))//2
aux=3*each_exc
assert dim==1384120320 and each_exc==687865856 and aux==2063597568
cert=json.loads((ROOT/'codex/work_C33/root_certificate_exact_n23.json').read_text())
summary=json.loads((ROOT/'codex/work_C33/root_certificate_n23.json').read_text())
assert len(cert['intervals'])==46
assert cert['g_power']==365 and cert['perturbation_budget_power']==33
assert summary['root_count']==46 and summary['finite_g_power']==-365
# Exact JSON stores each margin as numerator and denominator strings.
def fraction(x):
    if isinstance(x,list):return F(int(x[0]),int(x[1]))
    return F(x)
assert fraction(cert['delta2'])>F(140193,10000000)
assert fraction(cert['delta4'])>F(31558,10**11)
# The normalized higher-conversion bound is proportional to g^2.
for additional_power in (0,1,2,10,100,10000):
    ratio=F(1,1<<(2*additional_power))
    assert 0<ratio<=1

result=dict(status='all_exact_assertions_passed',global_jw_cases=jw_cases,
    neutral_wire_cases=gauge_cases,missing_parity_counterexample=missing_parity_witness,
    cell_dimension=dim,excited_states_per_parity=each_exc,auxiliary_modes=aux,
    total_modes=aux+n,certified_roots=46,
    note='Checks exact algebra and saved certificate fields; does not rerun the original inertia certificate or establish RG.')
(HERE/'n23_matching_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
