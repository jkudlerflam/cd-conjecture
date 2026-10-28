"""Calibrate physical-basis entropy against an exact one-component free sea."""
import json
import numpy as np
from xx_vacancy import HERE,build,ground,entropy
out=[]
for L in [8,12,16]:
    N=L//2; codes,mats=build(L,N,0)
    ev,vs,res=ground(-mats[0]); psi=vs[:,0]
    ks=2*np.pi/L*(np.arange(N)-(N-1)/2)
    errs=[]; ents={}
    for ell in range(1,L//2+1):
        exact=np.exp(1j*np.subtract.outer(np.arange(ell),np.arange(ell))[:,:,None]*ks).sum(axis=2)/L
        eig=np.linalg.eigvalsh(exact).clip(1e-15,1-1e-15)
        sexact=float(-np.dot(eig,np.log(eig))-np.dot(1-eig,np.log(1-eig)))
        sed=entropy(codes,psi,L,ell,N,0)
        errs.append(abs(sed-sexact));ents[str(ell)]=sed
    item={'L':L,'dimension':len(codes),'E0':float(ev[0]),'energy_error':abs(float(ev[0])+2/np.sin(np.pi/L)),
          'max_entropy_error':max(errs),'residual':res[0],'entropy':ents}
    assert item['energy_error']<1e-10 and item['max_entropy_error']<1e-9
    out.append(item)
coefficient=3*np.polyfit(np.log([d['L'] for d in out]),[d['entropy'][str(d['L']//2)] for d in out],1)[0]
result={'status':'numerical_independent_free_charge_calibration','cases':out,'half_chain_size_c':float(coefficient)}
(HERE/'calibration.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
