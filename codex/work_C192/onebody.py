# One-body exponents from the saved window matrices (run from codex/work_C176).
# g(r) = G(c - r//2, c - r//2 + r) for pairs centred on the chain centre (averaged over nearby centres).
# Demodulate one valley, low-pass with a binomial filter, fit |a(r)| ~ tan(pi r / 2L)^(-eta).
import json,numpy as np,sys
def last(f): return json.loads(open(f).read().strip().split('\n')[-1])
def eta(G,L,k,rlo,rhi,nc=6):
    W=G.shape[0]; c0=W//2; rs=np.arange(1,W-2*nc-2); g=[]
    for r in rs:
        g.append(np.mean([G[c-r//2,c-r//2+r] for c in range(c0-nc,c0+nc+1) if c-r//2>=0 and c-r//2+r<W]))
    h=np.array(g)*np.exp(-1j*k*rs)
    for _ in range(2): h=np.convolve(h,[.25,.5,.25],mode='same')
    s=(rs>=rlo)&(rs<=rhi); x=np.log(np.tan(np.pi*rs[s]/(2*L))); y=np.log(np.abs(h[s]))
    p=np.polyfit(x,y,1); return -p[0], np.std(y-np.polyval(p,x))
k1=float(sys.argv[1])
for tag in sys.argv[2:]:
    rec=last(f'typhon_out/out/{tag}.jsonl'); m=rec['meas']; L=rec['cfg']['L']
    out=[f'{tag:13s}']
    for nm in ('b','f'):
        G=np.array(m[f'G{nm}_re'])+1j*np.array(m[f'G{nm}_im']); W=G.shape[0]
        res=[eta(G,L,k1,lo,hi) for (lo,hi) in ((6,W//4),(10,W//3),(16,W//2),(24,int(.7*W)))]
        out.append(f'eta_{nm}: '+' '.join(f'{e:.3f}({s:.2f})' for e,s in res))
    print('  '.join(out))
