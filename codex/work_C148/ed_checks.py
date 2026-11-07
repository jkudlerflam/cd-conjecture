"""ED checks for Delta>1: (1) Bethe = ED sector minima at flux pi per bond, L<=12;
(2) end-to-end L(E_2L-2E_L)/M with field and DM, L<=6; (3) half-filling parity and h_c1."""
import json, itertools, numpy as np
from check_c129 import energy, bond_M
from circle import Circ, solve, Restricted
deltas=[1.001,1.01,1.1,1.5,2.,3.,5.,10.,20.]
res={}
err=0.; worst=None
for d in deltas:
    for L in range(2,13):
        for n in range(1,L//2+1):
            e=abs(float(solve(L,n,d)[0])-energy(L,n,d))
            if e>err: err,worst=e,(d,L,n)
res['bethe_vs_ed']=dict(max_error=err,at=worst,cases=len(deltas)*sum(L//2 for L in range(2,13)))
print(res['bethe_vs_ed'],flush=True)
# (2) end to end with field and DM: E_L(B) = min_n e_L(n, L theta) - B(L/2-n), theta=atan(D/Jperp)
rows=[]; mx=(-1,None)
for d in deltas:
    for dm in [0.,0.3,1.0]:
        J=1.0; Jperp=1/np.hypot(1,dm); D=dm*Jperp; theta=np.arctan2(D,Jperp); jz=d*J
        M0=Circ(d)
        for L in range(2,7):
            small={n:energy(L,n,d,theta,J) for n in range(L+1)}
            large={n:energy(2*L,n,d,2*theta/2,J) if False else energy(2*L,n,d,theta,J) for n in range(2*L+1)}
            hc2=J*(1+d)
            for B in np.r_[np.linspace(0,hc2*1.05,43)]:
                EL=min(small[n]-B*(L/2-n) for n in small)
                E2=min(large[n]-B*(L-n) for n in large)
                M=bond_M(J,jz,B)
                Qv=L*(E2-2*EL)/M
                rows.append(Qv)
                if Qv>mx[0]: mx=(Qv,dict(delta=d,dm=dm,L=L,B=float(B)))
res['end_to_end']=dict(samples=len(rows),max_ratio=mx[0],at=mx[1],min_ratio=float(min(rows)))
print(res['end_to_end'],flush=True)
# (3) parity at half filling and h_c1
par=[]
for d in [1.01,1.1,1.5,2.,5.]:
    M=Circ(d)
    for L in [9,10,33,34,129,130]:
        n=L//2
        E,_=solve(L,n,d); E2,_=solve(2*L,2*n,d)
        par.append(dict(delta=d,L=L,n=n,defect_over_M0=float(L*(E2-2*E)/M0.M0) if False else float(L*(E2-2*E)/M.M0)))
    # h_c1 = -e_bulk'(1/2-) = -h(pi/2) = 2 a r0(pi/2)
    hc1=2*M.a*float(M.r0(np.array([np.pi/2]))[0])
    gaps=[float(solve(L,L//2-1,d)[0]-solve(L,L//2,d)[0]) for L in [64,256,1024]]
    par.append(dict(delta=d,hc1_formula=hc1,finite_gaps_L64_256_1024=gaps))
res['parity_hc1']=par
for p in par: print(p)
json.dump(res,open('ed_checks.json','w'),indent=1)
