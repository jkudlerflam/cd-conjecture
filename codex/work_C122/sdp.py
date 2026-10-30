"""Conditional-min-entropy SDP, real charge-invariant matrices, NumPy/SciPy only.

Primal log-det path with Newton steps. Every returned interval uses independently
repaired primal and dual feasible points, not optimizer status alone.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigvalsh


def charges(d, n):
    digits = (np.arange(d**n)[:, None] // d**np.arange(n)) % d
    return [tuple(np.count_nonzero(row == a) for a in range(d-1)) for row in digits]


class MinSDP:
    def __init__(self, d, n, seam):
        self.d, self.n, self.seam = d, n, seam
        self.q, self.m = d**len(seam), d**(n-len(seam))
        rest = [i for i in range(n) if i not in seam]
        order = rest + seam
        digits = (np.arange(d**n)[:, None] // d**np.arange(n)) % d
        self.perm = digits @ (d**np.array(order))
        cq, dq = charges(d, len(seam)), charges(d, len(rest))
        basis = []
        for i in range(self.m):
            for j in range(i, self.m):
                if dq[i] == dq[j]:
                    b = np.zeros((self.m, self.m))
                    b[i,j] = b[j,i] = 1.
                    basis.append(b)
        self.basis = np.array(basis)
        self.cost = np.trace(self.basis, axis1=1, axis2=2)
        labels = [tuple(a+b for a,b in zip(c,x)) for c in cq for x in dq]
        self.blocks = []
        for label in sorted(set(labels)):
            ix = np.array([i for i,x in enumerate(labels) if x == label])
            ci, di = ix//self.m, ix%self.m
            f = self.basis[:,di[:,None],di[None,:]]*(ci[:,None] == ci[None,:])
            self.blocks.append((ix, f))

    def solve(self, matrix, tolerance=2e-6):
        w = matrix[np.ix_(self.perm, self.perm)]
        wb = [w[np.ix_(ix,ix)] for ix,_ in self.blocks]
        scale = max(float(eigvalsh(a)[-1]) for a in wb)
        x = self.cost*(scale+1.)
        dim = len(w)
        iterations = 0
        def evaluate(z, t, derivatives=True):
            val = float(self.cost@z)
            grad = self.cost.copy()
            hess = np.zeros((len(z),len(z)))
            inverses = []
            try:
                for (_, f), a in zip(self.blocks, wb):
                    slack = np.einsum('k,kij->ij',z,f)-a
                    cf = cho_factor(slack,lower=True,check_finite=False)
                    val -= t*2*np.log(np.diag(cf[0])).sum()
                    if derivatives:
                        inv = cho_solve(cf,np.eye(len(a)),check_finite=False)
                        inverses.append(inv)
                        grad -= t*np.einsum('ij,kji->k',inv,f)
                        prod = inv@f
                        hess += t*np.einsum('kij,lji->kl',prod,prod)
            except np.linalg.LinAlgError:
                return None
            return val, grad, hess, inverses
        for t in np.geomspace(.1, tolerance/(10*dim), 12):
            for _ in range(90):
                iterations += 1
                val,g,h,ivs = evaluate(x,t)
                step = np.linalg.solve(h,-g)
                descent = float(g@step)
                if -descent < max(1e-16,t*1e-9):
                    break
                rate = 1.
                for _ in range(60):
                    trial = evaluate(x+rate*step,t,False)
                    if trial is not None and trial[0] <= val+.01*rate*descent:
                        break
                    rate *= .5
                x += rate*step
            _,g,h,ivs = evaluate(x,t)
            omega = np.einsum('k,kij->ij',x,self.basis)
            mins = [float(eigvalsh(np.einsum('k,kij->ij',x,f)-a)[0])
                    for (_,f),a in zip(self.blocks,wb)]
            repair = max(0.,-min(mins))+1e-12*max(1.,scale)
            upper = float(np.trace(omega)+self.m*repair)
            dual = np.zeros_like(w)
            for (ix,_),inv in zip(self.blocks,ivs):
                dual[np.ix_(ix,ix)] = t*inv
            dual = (dual+dual.T)/2
            # Positive by construction, but repair any floating-point residue.
            mineig = min(float(eigvalsh(dual[np.ix_(ix,ix)])[0]) for ix,_ in self.blocks)
            dual += max(0.,-mineig)*np.eye(dim)
            trc = np.trace(dual.reshape(self.q,self.m,self.q,self.m),axis1=0,axis2=2)
            cap = max(1.,float(eigvalsh(trc)[-1]))*(1+1e-11)
            dual /= cap
            lower = float(np.sum(dual*w.T))
            if upper-lower < tolerance*max(1.,lower):
                break
        assert lower > 0 and upper >= lower-1e-9
        assert upper-lower < tolerance*max(1.,lower)*3, (lower,upper,np.linalg.norm(g))
        return dict(lower=lower,upper=upper,value=(lower+upper)/2,
                    gap=upper-lower,iterations=iterations,
                    primal_min_slack=min(mins)+repair,
                    dual_min_eigenvalue=mineig,
                    dual_partial_trace_max=float(eigvalsh(trc/cap)[-1]))
