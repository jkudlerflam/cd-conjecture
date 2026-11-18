"""C200: half filling, trial y^+ = quantiles of Z_0 at levels I_j/L + kappa0/L with kappa0 possibly > 1/2;
nodes whose level leaves the range of Z_0 are placed at +infinity.  Evaluates the Bethe map directly.
Prints min_j F_j/(2 pi eta) over finite nodes, and min over infinite nodes.  Floating point only."""
import numpy as np, sys
def run(eta, L, k0):
    n = L//2; g = np.pi-eta; b = np.pi/g
    j = np.arange(1, n+1); I = j-(n+1)/2
    lev = I/L+k0/L                      # Z_0 level, range (-1/4, 1/4)
    fin = lev < 0.25
    y = (2/b)*np.arctanh(np.tan(np.pi*lev[fin]))
    nf = fin.sum(); ninf = n-nf
    th2 = np.tan(eta/2); ct = 1/np.tan(eta)
    p = 2*np.arctan(np.tanh(y)*th2)
    F = L*p-2*np.pi*I[fin]
    # theta(x) = -2 atan(tanh(x) cot(eta))
    for s in range(0, nf, 512):
        d = y[s:s+512, None]-y[None, :]
        F[s:s+512] -= np.sum(-2*np.arctan(np.tanh(d)*ct), axis=1)
    F -= ninf*(np.pi-2*eta)             # theta(-inf) = pi - 2 eta for each node at +inf
    # infinite nodes, i-th from top, ordered among themselves
    Finf = np.array([eta*(L-2*(n-2*i+1)) for i in range(1, ninf+1)]) if ninf else np.array([np.inf])
    jm = np.argmin(F)
    return F.min()/(2*np.pi*eta), jm+1, nf, Finf.min()/(2*np.pi*eta)
if __name__ == '__main__':
    for eps in [0.5, 1e-1, 1e-2, 1e-4, 1e-6, 1e-8]:
        eta = np.arccos(1-eps)
        for L in [4, 16, 64, 256, 1024, 4096]:
            out = []
            for k0 in [0.4, 1.0, 2.0, 5.0, 20.0]:
                if k0 >= L/2: out.append('   -   '); continue
                m, jm, nf, mi = run(eta, L, k0)
                out.append(f'{m:8.4f}@{jm}')
            print(f'1+D={eps:7.1e} L={L:5d} ', '  '.join(out), flush=True)
