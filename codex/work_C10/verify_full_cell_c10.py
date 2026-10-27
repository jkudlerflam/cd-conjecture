"""Exhaustive quantum block ED for every computation string in the four-site cell.

Boundary computation markers and their data are fixed correctly. The selector
walker is factored out using the exact commutation proved in RESULT.md.
"""
import itertools
import json
import numpy as np
from scipy.linalg import eigvalsh, block_diag
from build_bond_c10 import make_bond, COMP, ALLOWED

J,kappa,Delta,g=100.,20.,1.,.2
C=3*abs(g)
I2=np.eye(2,dtype=complex)
I16=np.eye(16,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
omega=np.array([1,0,0,1],complex)/np.sqrt(2)
PO=np.outer(omega,omega.conj())
SW=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)

def embed(op,site):
    return np.kron(np.eye(2**site),np.kron(op,np.eye(2**(2-site))))

rules=[(('H0','R3'),('A2','H1'),SW),(('A2','H1'),('H2','R3'),np.eye(4))]
energies=[(('A1','H0'),3*Delta*(np.eye(4)-PO)),
          (('H0','R3'),1.5*g*(np.kron(X,X)-np.kron(Y,Y))),
          (('A1','H2'),1.5*g*(np.kron(X,X)+np.kron(Y,Y)))]
embedded_rules=[[embed(rule[2],site) for site in range(3)] for rule in rules]
embedded_energies=[[embed(term[1],site) for site in range(3)] for term in energies]
words=list(itertools.product(COMP,repeat=4))
wordid={word:i for i,word in enumerate(words)}
adj=[[] for _ in words]
violations=[]
for wi,word in enumerate(words):
    seq=('M',)+word+('M',)
    violations.append(sum((seq[j],seq[j+1]) not in ALLOWED for j in range(5))+word.count('M'))
    for site in range(3):
        pair=word[site:site+2]
        for ri,(a,b,W) in enumerate(rules):
            if pair==a:
                new=b
            elif pair==b:
                new=a
            else:
                continue
            dest=word[:site]+new+word[site+2:]
            adj[wi].append((wordid[dest],ri,site))

seen=set()
component_count=0
largest=0
illegal_min=np.inf
legal_spectrum=None
legal_comp=None
minimum_bad_word=None
component_sizes={}
for root in range(len(words)):
    if root in seen: continue
    pending=[root]
    seen.add(root)
    component=[]
    while pending:
        wi=pending.pop()
        component.append(wi)
        for wj,_,_ in adj[wi]:
            if wj not in seen:
                seen.add(wj)
                pending.append(wj)
    component_count+=1
    largest=max(largest,len(component))
    component_sizes[len(component)]=component_sizes.get(len(component),0)+1
    at={wi:j for j,wi in enumerate(component)}
    mat=np.zeros((16*len(component),16*len(component)),complex)
    is_legal=all(violations[wi]==0 for wi in component)
    assert is_legal or all(violations[wi]>0 for wi in component)
    for wi in component:
        word=words[wi]
        a=16*at[wi]
        diagonal=(kappa*violations[wi]+C*sum(x.startswith('H') for x in word)+(J/2)*len(adj[wi]))*I16
        diagonal=diagonal.copy()
        for site in range(3):
            pair=word[site:site+2]
            for ei,(pattern,energy) in enumerate(energies):
                if pair==pattern:
                    diagonal+=embedded_energies[ei][site]
        mat[a:a+16,a:a+16]=diagonal
        for wj,ri,site in adj[wi]:
            b=16*at[wj]
            mat[a:a+16,b:b+16]+=-(J/2)*embedded_rules[ri][site]
    ev=eigvalsh(mat)
    if is_legal:
        assert legal_spectrum is None
        legal_spectrum=ev
        legal_comp=component
    elif ev[0]<illegal_min:
        illegal_min=float(ev[0])
        minimum_bad_word=words[root]
assert component_count>0 and len(seen)==12**4
assert illegal_min>=kappa-1e-9

# Independently compress the actual 9216-by-9216 bond to the legal cell basis.
h,hs,hc,metadata=make_bond(J=J,kappa=kappa,Delta=Delta,g=g)
legal_clock=[('A1','H0','R3','R4'),('A1','A2','H1','R4'),('A1','H2','R3','R4')]
# Choose a nonmarker selector letter throughout the interior. hc is blind to its value.
sitebasis=[]
for word in legal_clock:
    for data in itertools.product(range(2),repeat=4):
        interior=[(1*12+COMP.index(c))*2+d for c,d in zip(word,data)]
        sitebasis.append(tuple([0]+interior+[0]))
actual=np.zeros((48,48),complex)
for ai,a in enumerate(sitebasis):
    for bi,b in enumerate(sitebasis):
        for bond in range(5):
            if any(a[s]!=b[s] for s in range(6) if s not in (bond,bond+1)):
                continue
            actual[ai,bi]+=hc[a[bond]*96+a[bond+1],b[bond]*96+b[bond+1]]
Ws=embed(SW,1)
Ud=block_diag(I16,Ws,Ws)
lap=np.array([[1,-1,0],[-1,2,-1],[0,-1,1]],float)/2
D0=embedded_energies[0][0]+embedded_energies[1][1]
D2=embedded_energies[2][0]
expected=J*Ud@np.kron(lap,I16)@Ud.conj().T+block_diag(D0,np.zeros((16,16)),D2)+(C-metadata['rho'])*np.eye(48)
projection_error=float(np.linalg.norm(actual-expected,2))
assert projection_error<1e-11
report=dict(interior_track_states=len(words),quantum_dimension=16*len(words),
            connected_components=component_count,component_size_histogram=component_sizes,
            largest_quantum_block=16*largest,
            legal_ground_energy_including_C=float(legal_spectrum[0]),
            legal_ground_dimension=int(np.sum(np.abs(legal_spectrum-legal_spectrum[0])<1e-8)),
            minimum_illegal_energy_including_C=illegal_min,
            minimum_illegal_energy_after_counterterm=illegal_min-metadata['rho'],
            analytic_illegal_bound_fixed_marker_labels=kappa-metadata['rho'],
            legal_ground_energy_after_counterterm=float(legal_spectrum[0]-metadata['rho']),
            an_illegal_minimizing_word=list(minimum_bad_word),
            actual_bond_projection_error=projection_error,
            scope='Every quantum block was diagonalized for all 12^4 interior computation strings with correct fixed endpoint marker labels and marker data. Endpoint-marker mismatches and arbitrary cell lengths are covered analytically, not by this enumeration.')
with open('full_cell_verification.json','w') as f:
    json.dump(report,f,indent=2)
print(json.dumps(report,indent=2))
