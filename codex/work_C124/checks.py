"""C124 checks. Python/numpy/scipy only. Outputs stay beside this file."""
from pathlib import Path
from itertools import product
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp
import xx_kernels as xx
import ed_kernels as ed

ROOT = Path(__file__).resolve().parent


def onsite(L, a, N):
    eps = L**-2
    tau = 2*a*L
    z = lambda t: float(logsumexp(-t*np.array([0., eps, 1.])))
    C = z(2*tau)+(2*L-2)*z(tau)
    G = z(2*N*tau)-2*z(N*tau)
    assert -np.log(3)-1e-12 <= G <= 1e-12
    return dict(L=L, a=a, N=N, epsilon=eps, boundary=C,
                boundary_per_L=C/L, G=G, rate=0.)


def histories_dense(h, d, n, beta, N, sites):
    ev, vec = eigh(h)
    f = (vec*np.exp(-beta*ev/(4*N)))@vec.conj().T
    e = f@f
    labels = (np.arange(d**n)[:, None]//(d**np.arange(n))[None, :])%d
    vals = []
    for history in product(range(d), repeat=N*len(sites)):
        history = np.array(history).reshape(N, len(sites))
        w = f.copy()
        for j in range(N):
            keep = np.all(labels[:, sites] == history[j], axis=1)
            w = keep[:, None]*w
            w = (f if j == N-1 else e)@w
        vals.append(float(np.sum(np.abs(w)**2)))
    return np.array(vals)


def embed_bond(k, n, d, i, j):
    states = np.arange(d**n)
    x = (states//d**i)%d
    y = (states//d**j)%d
    out = np.zeros((d**n, d**n), complex)
    for a, b in product(range(d), repeat=2):
        changed = states+(a-x)*d**i+(b-y)*d**j
        out[changed, states] += k[a*d+b, x*d+y]
    return out


def commuting_checks():
    # Each physical site is (left qubit, right qubit). The bond acts
    # only on right_i and left_(i+1). Adjacent bonds have disjoint support.
    x = np.array([[0,1],[1,0]], complex)
    y = np.array([[0,-1j],[1j,0]], complex)
    z = np.diag([1.,-1.])
    k0 = -.4*np.kron(x,x)-.23*np.kron(z,z)-.11*np.kron(z,np.eye(2))
    k1 = k0+.19*np.kron(y,x)+.07*np.kron(x,z)
    rows = []
    for name, k in [('stoquastic', k0), ('complex', k1)]:
        # Build the d=4 bond in order (left_i,right_i,left_j,right_j).
        h = np.zeros((16,16), complex)
        for li,ri,lj,rj,ro,lo in product(range(2), repeat=6):
            h[(2*li+ro)*4+2*lo+rj,(2*li+ri)*4+2*lj+rj] = k[2*ro+lo,2*ri+lj]
        comm = np.kron(h,np.eye(4))@np.kron(np.eye(4),h)
        comm -= np.kron(np.eye(4),h)@np.kron(h,np.eye(4))
        assert np.linalg.norm(comm)<1e-12
        H = {n:sum((embed_bond(h,n,4,i,(i+1)%n) for i in range(n)),
                   np.zeros((4**n,4**n),complex)) for n in [2,4]}
        for beta in [.2, 2., 6.]:
            for N in [1,2]:
                q = histories_dense(H[2],4,2,beta,N,[0])
                qd = histories_dense(H[4],4,4,beta,N,[0,2])
                # Explicit positive separator built from one-edge histories.
                fR = histories_dense(k,2,2,beta,N,[1])
                fL = histories_dense(k,2,2,beta,N,[0])
                histories = list(product(range(4),repeat=N))
                B = np.empty((4**N,4**N))
                for si,s in enumerate(histories):
                    ri = sum((s[j]%2)*2**(N-1-j) for j in range(N))
                    for ti,t in enumerate(histories):
                        li = sum((t[j]//2)*2**(N-1-j) for j in range(N))
                        B[si,ti] = fR[ri]*fL[li]
                # qd enumeration is time-major, with two sites per time.
                predicted = []
                for labels in product(range(4),repeat=2*N):
                    s = sum(labels[2*j]*4**(N-1-j) for j in range(N))
                    t = sum(labels[2*j+1]*4**(N-1-j) for j in range(N))
                    predicted.append(B[s,t]*B[t,s])
                err = max(np.max(abs(q-np.diag(B)))/max(q),
                          np.max(abs(qd-predicted))/max(qd))
                G = float(np.log(q@q)-np.log(qd.sum()))
                R = float(2*logsumexp(-beta*eigh(H[2],eigvals_only=True))
                          -logsumexp(-beta*eigh(H[4],eigvals_only=True)))
                assert err<1e-10,(name,beta,N,err)
                assert G<1e-10 and abs(R)<1e-10
                rows.append(dict(model=name,beta=beta,N=N,G=G,R=R,error=float(err)))
    return rows


def main():
    start=time.time()
    out = dict(commuting=commuting_checks(),
               onsite=[onsite(L,a,N) for L in [4,16,64,256,1024,8192]
                       for a in [.25,1.,4.] for N in [1,2,4,8]])
    print('commuting and onsite checks done', time.time()-start, flush=True)
    # Calibration of copied finite-grid formulas against full Hilbert matrices.
    out['calibration']=[]
    for n in [4,6]:
        ev,vec=xx.dense_xx(n)
        for beta in [.4,4.]:
            for N in [1,2]:
                for sites in [[0],[0,n//2]]:
                    q=xx.dense_histories(ev,vec,beta,N,sites)
                    qf=xx.xx_histories(n,beta,N,sites)
                    err=float(np.max(abs(q-qf))/max(q))
                    assert err<2e-10
                    out['calibration'].append(dict(n=n,beta=beta,N=N,sites=sites,error=err))
    out['xx']=[]
    for L in [4,16,128,1024,8192]:
        for a in [.25,1.,4.]:
            for N in [1,2,4]:
                row=xx.xx_case(L,a,N)
                # Inherited proved bound; no interacting extrapolation.
                upper=2*np.log(2)+(2*np.log(2)+(10*np.pi/3)*a)*N
                assert row['G']<=upper+1e-9
                out['xx'].append(row)
        print('XX',L,time.time()-start,flush=True)
    out['xx'].append(xx.xx_case(8192,1.,6))
    out['ed']=[]
    for model,d,delta,M,lengths in [
        ('XX',2,0.,.5,[4,6]),('XXZ_m09',2,-.9,.5,[4,6]),
        ('XXZ_p05',2,.5,.5,[4,6]),('Heisenberg',2,1.,.5,[4,6]),
        ('SU3',3,1.,1.,[3])]:
        for L in lengths:
            small,big=ed.sectors(L,d,delta),ed.sectors(2*L,d,delta)
            for a,N in [(a,N) for a in [.25,1.,4.] for N in [1,2,4]]+[(1.,8)]:
                row=ed.one_case(small,big,a,N,M,model)
                if model=='XX' and N<=4:
                    row['gaussian_error']=abs(row['G']-xx.xx_case(L,a,N)['G'])
                    assert row['gaussian_error']<2e-9
                out['ed'].append(row)
            print(model,L,'max block',big['maxblock'],time.time()-start,flush=True)
            out['seconds']=time.time()-start
            (ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    out['seconds']=time.time()-start
    (ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('complete',out['seconds'],flush=True)


if __name__=='__main__':
    main()
