"""Distance dependence of the pair correlator at fixed momenta (pi and the incommensurate peaks), run from meas/ (2026-11-24, Claude Opus 5.5)."""
import json, gzip, sys, numpy as np
def comp(v, r, q, w):
    out=[]
    for rc in [16,24,32,48,64,96,128,160]:
        s=(r>=rc-w)&(r<=rc+w)
        if s.sum()<2*w: continue
        out.append((rc, abs(np.mean(v[s]*np.exp(-1j*q*r[s])))))
    return np.array(out)
for fn,qs in [("T1_L384.json.gz",[1.0,0.793,1.207,0.207]),("T1_L256.json.gz",[1.0,0.793,1.207,0.207]),("T2_L384.json.gz",[1.0,0.882,1.118]),("T3_L384.json.gz",[1.0,0.793,1.207,0.207])]:
    d=json.load(gzip.open(fn))['meas']; L=d['L']
    print("==",fn)
    for d_ in (2,3):
        k="%d_%d"%(d_,L//2-3); v=np.array(d['Dpair'][k]); v=(v[:,0]+1j*v[:,1]); r=np.arange(len(v))+d_+1; v=v*(-1.0)**r
        for q in qs:
            e=comp(v,r,q*np.pi,10)
            sl=-np.polyfit(np.log(e[:,0]),np.log(e[:,1]),1)[0] if len(e)>2 else np.nan
            m=e[:,0]<=L/4
            sl2=-np.polyfit(np.log(e[m,0]),np.log(e[m,1]),1)[0]
            print("  d=%d q=%.3fpi  amp(r=16..160): %s  slope(all)=%.2f slope(r<=L/4)=%.2f"%(d_,q," ".join("%.1e"%x for x in e[:,1]),sl,sl2))
