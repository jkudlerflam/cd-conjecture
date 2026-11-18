"""C200: half filling, true Bethe roots against Theorem H.
(1) trial kappa0 = 1/2 (N(y_j) = j/L, top node at infinity): min_j Ftilde_j/(2 pi eta) >= 0
(2) roots: x_j <= y_j, i.e. L N(x_j) <= j;  by symmetry L N(x_j) >= j - 1;  sup_v |L N_x - L N| <= 1
(3) eta L sup |sigma_L - rho| <= (2/pi) eta cot eta, and min sigma_L/rho on {lambda >= 2}
Floating point only."""
import numpy as np, json
from bethe_ferro import Model, solve
from half_trial import run
rows = []
for eps in [0.29, 1e-1, 1e-2, 1e-4, 1e-6, 1e-8]:
    eta = np.arccos(1-eps); g = np.pi-eta; b = np.pi/g; m = Model(eta)
    for L in [2, 4, 8, 16, 64, 256, 1024, 2048]:
        n = L//2
        if L >= 4:
            tr = run(eta, L, 0.5)[0]
        else: tr = float('nan')
        G, x, info = solve(L, n, eta)
        N = lambda v: np.arctan(np.exp(b*v))/np.pi
        tN = L*N(x); j = np.arange(1, n+1)
        up = np.max(tN-j); lo = np.min(tN-(j-1))        # need up <= 0, lo >= 0
        # sup |L N_x - L N|: at root j the right value is j - tN_j, the left value j-1-tN_j
        supc = max(np.max(np.abs(j-tN)), np.max(np.abs(j-1-tN)))
        # density on a grid
        X = np.linspace(-min(x[-1]+3, 40), min(x[-1]+3, 40), 20001)
        T = lambda u: np.sin(2*eta)/(np.pi*(np.cosh(2*u)-np.cos(2*eta)))
        a1 = np.sin(eta)/(np.pi*(np.cosh(2*X)+np.cos(eta)))
        sig = a1.copy()
        for s in range(0, n, 256):
            sig += T(X[:, None]-x[None, s:s+256]).sum(axis=1)/L
        rho = b/(2*np.pi*np.cosh(b*X))
        dens = eta*L*np.max(np.abs(sig-rho)); bound = (2/np.pi)*eta/np.tan(eta)
        lam = eta*L*rho; sh = lam >= 2
        mr = float(np.min(sig[sh]/rho[sh])) if sh.any() else float('nan')
        rows.append(dict(eps=eps, L=L, trial_min=tr, up=float(up), lo=float(lo), supcount=float(supc), dens=float(dens), dens_bound=float(bound), min_ratio_shell=mr))
        print(f'1+D={eps:7.1e} L={L:5d} trial min F/(2pi eta)={tr:8.4f}  max(LN(x_j)-j)={up:8.4f}  min(LN(x_j)-(j-1))={lo:7.4f}  sup|LNx-LN|={supc:6.4f}  etaL|sig-rho|={dens:6.4f} <= {bound:6.4f}  min sig/rho on lam>=2: {mr:6.4f}', flush=True)
json.dump(rows, open('check_roots.json', 'w'), indent=1)
