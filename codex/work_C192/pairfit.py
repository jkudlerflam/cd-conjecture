import json,numpy as np,sys
from scipy.optimize import least_squares
def last(f): return json.loads(open(f).read().strip().split('\n')[-1])
k1=np.arcsin(0.25)
def fit(Dc,r,lo,hi,qs):
    s=(r>=lo)&(r<=hi); r=r[s]; y=Dc[s]; wt=r**0.8
    def model(al):
        B=np.stack([r**(-al[i])*np.exp(1j*q*r) for i,q in enumerate(qs)],1)
        A=np.linalg.lstsq(B*wt[:,None],y*wt,rcond=None)[0]
        return B,A
    def res(al):
        B,A=model(al); d=(B@A-y)*wt; return np.concatenate([d.real,d.imag])
    best=None
    for a0 in ([0.8]*len(qs),[1.5]*len(qs),[0.8,1.5,1.5][:len(qs)],[1.5,0.8,0.8][:len(qs)]):
        o=least_squares(res,a0,bounds=(0,6))
        if best is None or o.cost<best.cost: best=o
    B,A=model(best.x)
    return best.x,np.abs(A),np.sqrt(2*best.cost/np.sum(np.abs(y*wt)**2))
for tag in sys.argv[1:]:
    rec=last(f'typhon_out/out/{tag}.jsonl'); m=rec['meas']
    for key in sorted(m['Dpair']):
        d=int(key[0]); D=np.array(m['Dpair'][key]); Dc=D[:,0]+1j*D[:,1]; r=np.arange(len(Dc))+d+1.0
        for lo,hi in ((8,48),(12,72),(16,96)):
            al,A,rr=fit(Dc,r,lo,hi,[np.pi,2*k1,-2*k1])
            print(f'{tag} {key} r[{lo},{hi}] alpha_pi={al[0]:.2f} alpha_+2k1={al[1]:.2f} alpha_-2k1={al[2]:.2f} amps={A[0]:.3f},{A[1]:.3f},{A[2]:.3f} relres={rr:.3f}')
