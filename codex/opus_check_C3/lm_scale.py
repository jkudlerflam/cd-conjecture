import numpy as np, scipy.linalg as la
exec(open('lm_root.py').read().split('# Jacobian')[0].replace('for t in range(300)','for t in range(120)'))
hh=lambda x:-x*np.log(x)-(1-x)*np.log(1-x)
def expm_n(M):
    w,v=la.eigh((M+M.conj().T)/2); E=(v*np.exp(w-w.max()))@v.conj().T; Z=np.trace(E).real; return E/Z, np.log(Z)+w.max()
def kr(*ms):
    out=np.eye(1)
    for m in ms: out=np.kron(out,m)
    return out
I2=np.eye(2)
seen={}
for c,ev,_ in sols:
    key=round(c,5)
    if key in seen: continue
    seen[key]=True
# re-solve to get the sigma for each class
rng=np.random.default_rng(7); reps={}
for t in range(400):
    p0=rng.normal(size=16)*[0.3,1,3][t%3]
    r=least_squares(resid,p0,bounds=(-12,12),xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
    if np.abs(r.fun).max()<1e-10:
        s=sig(r.x); lam=step(s)[3]; k=round(3*lam/h14,5)
        if k not in reps: reps[k]=s
print("classes:",sorted(reps))
for k,s in sorted(reps.items()):
    r3,a,b,lam=step(s); t=ptr(s,[0],[2,2]); c=3*lam/h14
    L3=logm(r3); Ls=logm(s); Lt=logm(t)
    # (1,1,2): A=1,B=1,C=2 ; eta=1/3 ; (1-eta) log r4 = log s_AB + log r3_BCD - eta(log t_A + log s_CD) - (1-eta) log t_B + lam'
    eta=1/3
    X=kr(Ls,I2,I2)+kr(I2,L3)-eta*(kr(Lt,I2,I2,I2)+kr(I2,I2,Ls))-(1-eta)*kr(I2,Lt,I2,I2)
    r4,logZ=expm_n(X/(1-eta)); lam112=-(1-eta)*logZ
    m123=ptr(r4,[0,1,2],[2]*4); m234=ptr(r4,[1,2,3],[2]*4)
    err_ext=max(np.abs(m123-r3).max(),np.abs(m234-r3).max())
    print(f"class c_D={c:.5f}: lam(1,1,2)={lam112:.5f} vs (c/3)h(1/3)={c/3*hh(1/3):.5f}; 4-site marginal mismatch={err_ext:.3e}; min eig rho3={la.eigvalsh(r3).min():.3e}; offdiag(sigma)={np.abs(s-np.diag(np.diag(s))).max():.2e}")
