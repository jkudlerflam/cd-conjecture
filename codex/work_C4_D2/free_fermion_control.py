"""Exact finite-ring XX energies, including the parity-dependent boundary twist.

Independent verification by ED at L=8,10,12. The energy-fit failure is then
evaluated without an assumption about numerical eigensolver convergence.
"""
from d2 import quick, detailed, dump
import numpy as np

def ground(L,hz):
    energies=[]
    for N in range(L+1):
        k=2*np.pi*(np.arange(L)+(0.5 if N%2==0 else 0))/L
        eps=4*np.cos(k)-2*hz
        energies.append(hz*L+np.sort(eps)[:N].sum())
    N=int(np.argmin(energies))
    return dict(L=L,E0=float(energies[N]),N=N,m=0 if N%2==0 else L//2)

out={}
for hz in (.3,1.1):
    J=np.diag([1.,1.,0.]); field=np.array([0.,0.,hz]); rows=[]; errors=[]
    exact=[ground(L,hz) for L in range(8,23,2)]
    for r in exact[:3]: errors.append(abs(quick(J,field,r['L'])['E0']-r['E0']))
    assert max(errors)<1e-10
    L=np.array([r['L'] for r in exact],float); E=np.array([r['E0'] for r in exact])
    a=np.linalg.lstsq(np.column_stack([np.ones(len(L)),L**-2,L**-4]),E/L,rcond=None)[0][1]
    velocity=4*np.sqrt(1-(hz/2)**2)
    for r in exact:
        row=detailed(J,field,r['L'],m=r['m'],nup=r['N'])
        rows.append(row)
        print(hz,r['L'],row['c'],row['stress'][0]['selected']['c_MV'],flush=True)
    out[str(hz)]=dict(exact=exact,ED_energy_errors=errors,velocity=velocity,
                     naive_energy_c=float(-6*a/(np.pi*velocity)),rows=rows)
    dump('results/free_fermion_control.json',out)
