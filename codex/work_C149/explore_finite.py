import numpy as np, sys
from bethe_ferro import Model, residual, solve
from continuum import Kernel
from trial import continuum_trial, quantiles
eps=float(sys.argv[1]); eta=2*np.arcsin(np.sqrt(eps/2)); m=Model(eta); ker=Kernel(eta)
for L in [int(a) for a in sys.argv[2:]]:
  for n in [L//8, L//4, 3*L//8, L//2-2]:
    R,P=continuum_trial(eta,L,n)
    y=quantiles(P.inv,L,n)
    Phi=-residual(m,y,L)/(2*np.pi*L)
    rho=P.rho_at(y); 
    drho=np.interp(np.abs(y),R.y,np.gradient(R.r,R.y))*np.sign(y)
    lam=eta*L*rho
    pred=eta**2*drho/np.expm1(2*np.pi*lam)
    G,x,d=solve(L,n,eta)
    disp=L*(P.Z(x)-(np.arange(1,n+1)-(n+1)/2)/L)
    j=np.argmax(np.abs(Phi))
    print(f"L={L} n={n} Q={R.Q:.3f} lam0={eta*L*R.r[0]:.3g} lamQ={eta*L*R.r[-1]:.3g} max L|Phi|/eta={L*np.abs(Phi).max()/eta:.3f} at d={R.Q-abs(y[j]):.3g} lam={lam[j]:.3g}; max L|Phi-pred|/eta={L*np.abs(Phi-pred).max()/eta:.3f}; max|disp|={np.abs(disp).max():.3f}",flush=True)
