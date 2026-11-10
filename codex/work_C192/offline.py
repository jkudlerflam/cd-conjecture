# Two-valley exponents off the symmetric line (run from codex/, args: jsonl paths).
# Boson one-body function at each valley momentum, and the inter-valley pair function,
# with the open-chain conformal forms. Momenta are located from the data.
# M_aa = 2 eta_a, M_bb = 2 eta_b, M_ab = eta_pair - eta_a - eta_b, N = M^-1.
# Backscattering dimension N_aa + N_bb - 2 N_ab, two-pair vertex dimension 2 eta_pair,
# density structure factor slope * pi = sum(N)/2 (equals K_c on the line).
import json,numpy as np,sys
def last(f): return json.loads(open(f).read().strip().split('\n')[-1])
d=lambda x,L: np.sin(np.pi*x/(2*L))
def sm(h,n):
    for _ in range(n): h=np.convolve(h,[.25,.5,.25],mode='same')
    return h
def peaks(rs,g,npk,lo=-np.pi,hi=np.pi):
    ks=np.linspace(lo,hi,8001)[:-1]; w=np.cos(np.pi*(rs-rs[0])/(2*(rs[-1]-rs[0]+1)))**2
    F=np.abs(np.exp(-1j*np.outer(ks,rs))@(g*w)); out=[]
    for _ in range(npk):
        i=np.argmax(F); out.append((ks[i],F[i])); dk=np.abs(np.angle(np.exp(1j*(ks-ks[i])))); F=np.where(dk<0.6,0,F)
    return out
def fit(x,y):
    p=np.polyfit(x,y,1); return -p[0],np.std(y-np.polyval(p,x))
for f in sys.argv[1:]:
    rec=last(f); m=rec['meas']; L=rec['cfg']['L']; tag=f.split('/')[-1][:-6]
    G=np.array(m['Gb_re'])+1j*np.array(m['Gb_im']); W=G.shape[0]; c0=W//2; nc=6
    rs=np.arange(1,W-2*nc-2)
    g=np.array([np.mean([G[c-r//2,c-r//2+r] for c in range(c0-nc,c0+nc+1)]) for r in rs])
    xw=np.arange(W); ks=2*np.pi*np.arange(256)/256
    n2=np.array([(np.exp(-1j*k*xw)@G@np.exp(1j*k*xw)).real/W for k in ks]); w0,wpi=n2[np.cos(ks)>0].sum(),n2[np.cos(ks)<0].sum()
    pk=sorted(peaks(rs,g,2),key=lambda t:-(wpi if np.cos(t[0])<0 else w0))   # majority valley first
    frac=[(wpi if np.cos(k)<0 else w0)/(w0+wpi) for k,_ in pk]
    print(f'{tag}  L={L}  Wa,Wb={m["Wa_b"]:.2f},{m["Wb_b"]:.2f}  n(k) weight fractions maj,min={frac[0]:.3f},{frac[1]:.3f}')
    etas=[]
    kref=[]
    for lab,(k,amp) in zip(('maj','min'),pk):
        hh=sm(g*np.exp(-1j*k*rs),4); sl=(rs>=8)&(rs<=W//2); k=k+np.polyfit(rs[sl],np.unwrap(np.angle(hh[sl])),1)[0]; kref.append(k)
        h=np.abs(sm(g*np.exp(-1j*k*rs),4)); x=np.log(np.tan(np.pi*rs/(2*L))); y=np.log(h)
        res=[fit(x[(rs>=lo)&(rs<=hi)],y[(rs>=lo)&(rs<=hi)]) for lo,hi in ((6,W//4),(10,W//3),(16,W//2),(24,int(.7*W)))]
        etas.append(res); print(f'  boson {lab} k={k:+.4f} amp={amp:.2f}  eta: '+' '.join(f'{e:.3f}({s:.3f})' for e,s in res))
    pe={}
    for sep in ('2','3'):
        dat=[]
        for key,v in m['Dpair'].items():
            if not key.startswith(sep+'_'): continue
            x0=int(key.split('_')[1])+0.5+int(sep)/2; a=np.array(v); D=a[:,0]+1j*a[:,1]; r=np.arange(len(D))+int(sep)+1
            dat.append((x0,r,D))
        n0=min(len(D) for _,_,D in dat); Q=peaks(dat[0][1][:n0],np.mean([D[:n0] for _,_,D in dat],0),1,lo=np.pi-1,hi=np.pi+1)[0][0]
        sl=[np.polyfit(r[(r>=8)&(r<=int(.4*L))],np.unwrap(np.angle(sm(D*np.exp(-1j*Q*r),2)[(r>=8)&(r<=int(.4*L))])),1)[0] for _,r,D in dat]; Q=Q+np.mean(sl)
        res=[]
        for lo,hi in ((6,L//8),(10,L//5),(16,L//3),(24,int(.4*L))):
            X=[];Y=[]
            for x0,r,D in dat:
                h=np.abs(sm(D*np.exp(-1j*Q*r),2)); F=d(r,L)*d(2*x0+r,L)/np.sqrt(d(2*x0,L)*d(2*(x0+r),L))
                s=(r>=lo)&(r<=hi); X+=list(np.log(F[s])); Y+=list(np.log(h[s]))
            res.append(fit(np.array(X),np.array(Y)))
        pe[sep]=res; pe['Q'+sep]=Q; print(f'  pair d={sep} Q-pi={Q-np.pi:+.4f}  eta: '+' '.join(f'{e:.3f}({s:.3f})' for e,s in res))
    ksum=np.angle(np.exp(1j*(kref[0]+kref[1]-np.pi)))
    print(f'  momenta: k_maj+k_min-pi={ksum:+.4f}  Q-pi (d=2,3)={pe["Q2"]-np.pi:+.4f},{pe["Q3"]-np.pi:+.4f}')
    for wi in (2,3):
        ea,eb=etas[0][wi][0],etas[1][wi][0]; ep=0.5*(pe['2'][wi][0]+pe['3'][wi][0])
        M=np.array([[2*ea,ep-ea-eb],[ep-ea-eb,2*eb]]); N=np.linalg.inv(M)
        print(f'  window {wi+1}: eta_maj={ea:.3f} eta_min={eb:.3f} eta_pair={ep:.3f} M_ab={M[0,1]:+.3f}  bs_dim={N[0,0]+N[1,1]-2*N[0,1]:.2f}  O4_dim={2*ep:.2f}  sumN/2={N.sum()/2:.2f}')
