#!/usr/bin/env python3
"""Physical consistency checks on saved states' measured observables."""
import json
import numpy as np
from threadpoolctl import threadpool_limits
from common import HERE


def main():
    rows=[]
    with threadpool_limits(limits=1):
        for f in sorted((HERE/'outputs').glob('*.json')):
            d=json.loads(f.read_text())
            if d['status']!='done':continue
            L=d['cfg']['L'];P=d['cfg']['P'];m=d['meas']
            par_error=abs(d['parity']-(1-2*P));assert par_error<1e-8
            dens=np.array(m['dens']);assert min(dens)>-1e-9 and max(dens)<1+1e-9
            S=np.array(d['S']);ell=np.arange(1,L)
            assert min(S)>-1e-9 and max(S-np.minimum(ell,L-ell)*np.log(2))<1e-8
            G=np.array(m['Gf_re'])+1j*np.array(m['Gf_im']);ev=np.linalg.eigvalsh(G)
            herm=float(np.max(abs(G-G.conj().T)));assert herm<1e-12 and min(ev)>-1e-8 and max(ev)<1+1e-8
            lo,n=m['gwin'];diagerr=float(np.max(abs(G.diagonal().real-dens[lo:lo+n])));assert diagerr<1e-12
            rows.append(dict(tag=d['cfg']['tag'],parity_error=par_error,fermion_hermiticity_error=herm,
                fermion_eigenvalue_range=[float(min(ev)),float(max(ev))],fermion_diagonal_error=diagerr,
                density_reflection_error=float(np.max(abs(dens-dens[::-1])))))
    (HERE/'output_checks.json').write_text(json.dumps(dict(checked=len(rows),rows=rows,
        scope='Parity, entropy bounds and one-body positivity; not a CFT or ground-state certification.'),indent=2)+'\n')
    print('Passed physical-observable checks on',len(rows),'completed outputs.')


if __name__=='__main__':main()
