"""C164 calibration controls (Claude Opus agent, 2026-11-08).
dmrg  : open-chain DMRG c fits (dmrg_c.profile_fits) at a known c = 1 point (XX chain, a B_2 bond) and at a c = 2
        control outside B_2: two decoupled XX chains, written as one chain of d = 4 super-sites (a two-site cell of
        the qubit chain, i.e. couplings between qubits 2j, 2j+2 and 2j+1, 2j+3, which are range 3, not nearest
        neighbour). L counts sites of the respective chain (super-sites for the control).
fes   : two-site-cell damped VUMPS (vumps.py, C153/C158) for the same two models, output in the fesx.py format
        read by fesfit.py.
usage: controls.py dmrg MODEL L  |  controls.py fes MODEL chi1 chi2 ... out.json"""
import sys, json, time, os
import numpy as np
import est, dmrg_d as D
from dmrg_c import profile_fits

def h_xx():
    return est.bond_from_JB(np.diag([1.0, 1.0, 0.0]), np.zeros(3))

def h_two_xx():
    hx = h_xx().reshape(2, 2, 2, 2)                         # (s1 s2, t1 t2)
    one = np.eye(4).reshape(2, 2, 2, 2)
    # super-site = (a, b); bond acts on (a_j, b_j, a_{j+1}, b_{j+1})
    T = np.einsum('ACac,BDbd->ABCDabcd', hx, one) + np.einsum('ACac,BDbd->ABCDabcd', one, hx)
    return T.reshape(16, 16)

MODELS = dict(xx=h_xx, twoxx=h_two_xx)

if __name__ == '__main__':
    mode, model = sys.argv[1], sys.argv[2]
    h = MODELS[model]()
    if mode == 'dmrg':
        L = int(sys.argv[3]); t0 = time.time()
        chis = (16, 32, 64, 96, 128, 128) if L <= 64 else (16, 32, 64, 128, 192, 256, 256)
        r = D.dmrg(h, L, chis=chis)
        cp, cs, mid = profile_fits(r['S'], L)
        print(json.dumps(dict(model=model, L=L, c_prof=cp, c_sm=cs, Smid=mid, trunc=r['hist'][-1][2], t=time.time() - t0)))
    else:
        import vumps as V
        chis = [int(x) for x in sys.argv[3:-1]]; out = sys.argv[-1]
        d = int(round(np.sqrt(h.shape[0])))
        W = V.mpo_from_term(h, 2, d=d)
        data = dict(args=dict(n=2, model=model), results=[]); t0 = time.time()
        def on_chi(psi, rec):
            rec['time_s'] = time.time() - t0; data['results'].append(rec)
            json.dump(data, open(out + '.tmp', 'w'), default=float); os.replace(out + '.tmp', out)
        V.run([W] * 2, 2, chis, d=d, tol=1e-9, maxit=600, seed=0, mode='damped', on_chi=on_chi)
