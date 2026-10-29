"""Independent finite-matrix and factorized-state checks for C67."""
import ed_core as ed
from ed_z2 import couplings
import numpy as np
import json
from pathlib import Path

X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]])
Z=np.diag([1.,-1.]); I=np.eye(2)
P=[X,Y,Z]

def dense(L,J,b):
    def term(ops):
        out=np.ones((1,1))
        for j in reversed(range(L)): out=np.kron(out,ops.get(j,I))
        return out
    H=np.zeros((2**L,2**L),complex)
    for j in range(L):
        for a in range(3):
            H+=b[a]*term({j:P[a]})
            for c in range(3): H+=J[a,c]*term({j:P[a],(j+1)%L:P[c]})
    return H

out={}
v=[-.73,.42,.19,.31,-.27]; J,b=couplings(v); L=6
H=dense(L,J,b); full=np.linalg.eigvalsh(H); ee=[]; inter=[]
for p in [0,1]:
    data=ed.translation_orbits(L,p)
    for m in range(L):
        sec=ed.Sector(L,m,p,data); block=sec.matrix(J,b).toarray()
        W=np.column_stack([sec.expand(np.eye(len(block))[:,i]) for i in range(len(block))])
        inter.append(float(np.linalg.norm(H@W-W@block)))
        ee.extend(np.linalg.eigvalsh(block))
out['dense_spectrum_error']=float(np.max(np.abs(np.sort(ee)-full)))
out['intertwining_error']=max(inter)

A=1.; B=.5; C=.25; h=2*np.sqrt((A+C)*(B+C)); ct=np.sqrt((B+C)/(A+C))
theta=np.arccos(ct); spin=np.array([np.cos(theta/2),np.sin(theta/2)])
out['factorized']=[]
for D in [0.,.1,.2]:
    J,b=couplings([-A,-B,C,D,-h])
    H=dense(6,J,b); psi=spin
    for j in range(5): psi=np.kron(psi,spin)
    e=6*(-A-B-C); es=np.linalg.eigvalsh(H)
    out['factorized'].append(dict(D=D,h=h,cos_theta=ct,energy=e,
       residual=float(np.linalg.norm(H@psi-e*psi)),ground_error=float(es[0]-e),
       ground_multiplicity=int(np.sum(np.abs(es-e)<1e-9)),gap_above=float(es[2]-es[0])))

out['controls']=[]
for line in Path('controls.jsonl').read_text().splitlines():
    r=json.loads(line); L=r['L']
    if r['name']=='ising': exact=-2/np.sin(np.pi/(2*L))
    elif r['name']=='xx': exact=-4/np.sin(np.pi/L)
    else: continue
    out['controls'].append(dict(name=r['name'],L=L,error=abs(r['E0']-exact)))

Path('analytic_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
