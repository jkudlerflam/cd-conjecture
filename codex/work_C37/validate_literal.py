"""Small-size numerical audit against the supplied literal 23-letter cell.
This is a cross-check, not the large-n certificate algorithm.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,time
import numpy as np
sys.path.insert(0,'../../contrib/routeA_spatial')
from ycell import params,build_cell,code,second_order,logical
from structured import Leading,visibility
from modular import matrices,P
from certify_pairs import poly,pair_polynomial,gcd

start=time.monotonic();records=[]
for n in [1,2,3]:
    T0,M0,_=Leading(n,D=2,a=3,td=.03,loopflux=np.pi).calculate()
    r=[]
    for g in [.04,.02,.01]:
        pp=params(Dx=2,DM=2,Delta=3-2*np.cos(np.pi/(n+2)),delta=0,gp=g,gm=g,td=.03,
            phx={'0':np.pi/(2*n),'+':0.,'-':0.},phd={'0':np.pi,'+':0.,'-':0.})
        cell=build_cell(n,pp);E,G,gap=code(cell);H2,_=second_order(cell,G,E);h=logical(H2,n)
        r.append((h['T']/g**2,(h['muA']+h['muB'])/g**2))
    Told=(4*r[1][0]-r[0][0])/3;Mold=(4*r[1][1]-r[0][1])/3
    T=(4*r[2][0]-r[1][0])/3;M=(4*r[2][1]-r[1][1])/3
    sc=np.max(abs(T0))
    records.append(dict(n=n,cell_dimension=len(cell['H']),code_bond_dimension=len(H2),
      T_relative_error=float(np.max(abs(T-T0))/sc),
      M_error_over_T=float(np.max(abs(M-M0))/sc),
      M_imag_error_over_T=float(np.max(abs((M-M0).imag))/sc),
      M_real_error_over_T=float(np.max(abs((M-M0).real))/sc),
      T_previous_error=float(np.max(abs(Told-T0))/sc),
      M_previous_error=float(np.max(abs(Mold-M0))/sc),
      measured_gap=float(gap)))
controls=[]
for n in range(1,5):
    for label,tzero,ezero in [('tree',True,False),('zero_field',False,True)]:
        T,Tc,M,eps,ii=matrices(n,td_zero=tzero,eps_zero=ezero)
        f=poly(T,Tc,M,eps,P);q=pair_polynomial(f,P)
        gg=gcd(q,[(j*q[j])%P for j in range(1,len(q))],P)
        if n>=2 or ezero:assert len(gg)>1
        controls.append(dict(n=n,control=label,gcd_degree=len(gg)-1))
result=dict(literal=records,controls=controls,seconds=time.monotonic()-start)
json.dump(result,open('literal_results.json','w'),indent=2);print(json.dumps(result),flush=True)
