"""C120 small-ring utilities adapted from inspected C118 code."""
import numpy as np
from scipy.special import logsumexp
def xx_logz(n, beta):
    # Long double limits cancellation of the bulk term in large-ring ratios.
    answers = []
    for shift, parity in [(0.5, 0), (0.0, 1)]:
        eps = np.cos(2*np.longdouble(np.pi)*(np.arange(n, dtype=np.longdouble)+shift)/n)
        even, odd = np.longdouble(0), np.longdouble(-np.inf)
        for energy in eps:
            even, odd = np.logaddexp(even, odd-beta*energy), np.logaddexp(odd, even-beta*energy)
        answers.append((even, odd)[parity])
    return np.logaddexp(*answers)


def stats(ev, beta):
    shifted = ev-ev[0]
    lw = logsumexp(-beta*shifted)
    weights = np.exp(-beta*shifted-lw)
    mean = np.dot(weights, shifted)
    variance = np.dot(weights, (shifted-mean)**2)
    return float(lw-beta*ev[0]), float(lw), float(variance)


def local_bond(model):
    if model == 'SU3':
        return np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9), 3
    x = np.array([[0.,1.],[1.,0.]])/2
    y = np.array([[0.,-1j],[1j,0.]])/2
    z = np.diag([.5,-.5])
    h = np.kron(x,x)+np.kron(y,y)
    if model == 'Heisenberg':
        h += np.kron(z,z)
    return h.real, 2


def embed(h, d, n, edges):
    powers = d**np.arange(n, dtype=np.int64)
    result = np.zeros((d**n, d**n))
    for state in range(d**n):
        for i,j,weight in edges:
            a,b = (state//powers[i])%d, (state//powers[j])%d
            for aa in range(d):
                for bb in range(d):
                    value = weight*h[d*aa+bb,d*a+b]
                    if value:
                        target = state+(aa-a)*powers[i]+(bb-b)*powers[j]
                        result[target,state] += value
    return result


def entropy(rho):
    ev = np.linalg.eigvalsh((rho+rho.T)/2)
    ev = ev[ev>1e-15]
    return float(-np.dot(ev,np.log(ev)))


