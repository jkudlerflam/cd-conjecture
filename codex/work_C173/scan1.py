import numpy as np, json, sys
from eff_ed import *
out=open('scan1.jsonl','a')
for L in (24,28):
  for wa in (0.5,1.0,2.0):
    w=wa*np.exp(1j*1.3)
    for N in ([5,6,7,8] if L==24 else [6,7,8,9]):
      res={'L':L,'w_abs':wa,'phi':1.3,'N':N}
      for name,kw in (('blockade',{}),('blk_noV',{'V2':0,'V3':0}),('free',{'blockade':False,'V2':0,'V3':0})):
        if name=='free' and L>24: continue
        b,e,v=ground(L,N,w,k=1,**kw)
        S=entropies(b,v[:,0],L)
        c1,c2=cc_fit(S,L)
        res[name]={'E':float(e[0]),'c_fit':c1,'c_par':c2,'Smid':float(S[L//2-1])}
      print(json.dumps(res),flush=True); out.write(json.dumps(res)+'\n'); out.flush()
