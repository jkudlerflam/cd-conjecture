"""Dominant momenta of the pair, one-body and density correlators in meas/*.json.gz (2026-11-24, Claude Opus 5.5).
usage: python kpair.py ../meas/T1_L384.json.gz ..."""
import json, gzip, sys, numpy as np
def peaks(v, r, nq=8192, top=3):
    # v(r) complex; weight to suppress ends; return dominant q in [0,2pi) of sum v e^{-iqr}
    q = np.linspace(0, 2*np.pi, nq, endpoint=False)
    w = np.hanning(len(r))
    F = np.abs(np.exp(-1j*np.outer(q, r)) @ (v*w))
    out=[]; F2=F.copy()
    for _ in range(top):
        k=np.argmax(F2); out.append((q[k], F2[k]/F.max()))
        dq=2*np.pi/len(r)*2.5; F2[np.abs((q-q[k]+np.pi)%(2*np.pi)-np.pi)<dq]=0
    return out
for fn in sys.argv[1:]:
    d=json.load(gzip.open(fn))['meas']; L=d['L']
    print("==", fn, "L=%d"%L)
    for d_ in (2,3):
        for i0 in (L//2-3, L//2):
            k="%d_%d"%(d_,i0)
            if k not in d['Dpair']: continue
            v=np.array(d['Dpair'][k]); v=v[:,0]+1j*v[:,1]
            r=np.arange(len(v))+d_+1
            sel=(r>=L//16)&(r<=L//2)
            v_orig=v*(-1.0)**r          # ims frame -> original frame: pair phase (-1)^(j-i)
            pk=peaks(v_orig[sel], r[sel])
            print("  pair d=%d i0=%d: q/pi (orig frame) ="%(d_,i0), " ".join("%.4f(%.2f)"%(a/np.pi,b) for a,b in pk),
                  "  -> 2q/pi mod 2 = %.4f"%((2*pk[0][0]/np.pi)%2))
    # one-body, centred row i = L/2 - 1 ... j>i  (ims frame, report raw)
    i=L//2-40; re,im=d['Gb_rows'][str(i)]; g=np.array(re)+1j*np.array(im); r=np.arange(1,len(g)+1)
    sel=(r>=8)&(r<=L//2)
    print("  one-body row i=%d (ims frame raw): q/pi ="%i, " ".join("%.4f(%.2f)"%(a/np.pi,b) for a,b in peaks(g[sel],r[sel])))
    C=np.array(d['Cnn'][str(L//2)]); r=np.arange(1,len(C)+1); sel=(r>=8)&(r<=L//2)
    print("  density: q/pi =", " ".join("%.4f(%.2f)"%(a/np.pi,b) for a,b in peaks(C[sel].astype(complex),r[sel],top=4)))
