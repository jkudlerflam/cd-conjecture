"""(1) Compare the referee rebuild (my15.py) with C75's bond15.py entry by entry.
bond15.py is only imported here, as the object being checked, from a path
given on the command line (a scratch copy of codex/work_C75)."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.abspath(sys.argv[1]))
sys.dont_write_bytecode = True
import my15
import bond15, bond18

out = []
for n, kw in [(2, dict(K=100., g=.071, td=.003, flux=.34, theta=.23, field=.17, musel=.4, D=7.3, Delta=2.9)),
              (15, dict(K=2.**40, g=2**-8, td=1e-4, theta=-np.pi/4, field=1e-3, musel=1.9015866376241617, D=1000.))]:
    p = bond18.params(n, **kw)
    q = my15.default_params(n, D=p['D'], Delta=p['Delta'], g=p['g'], td=p['td'], phiA=p['flux'], phiB=0.,
                            psiA=n*p['flux']+np.pi/2, psiB=np.pi/2, portA=np.exp(1j*p['theta']),
                            portB=np.exp(-1j*p['theta']), K=p['K'], musel=p['musel'], field=p['field'])
    for lam in [1., 0.37]:
        H = bond15.bond(p, lam)
        # map bond15 letter tuples to referee names
        def name(l):
            k, s, a = l
            return {'L': f'L{a}', 'R': f'R{a}', 'd': f'D{s}{a}', 'x': f'X{s}{a}',
                    'e': 'M0', 'm': f'M{s}'}[k]
        perm = [my15.IX[name(l)] for l in bond15.LET]
        P = np.zeros((15, 15)); P[perm, range(15)] = 1
        PP = np.kron(P, P)
        Hc = PP @ H @ PP.T
        M = my15.bond(q, lam)
        q2 = np.array([my15.charge(l) for l in my15.LETTERS]); Q = (q2[:, None] + q2[None, :]).ravel()
        r = dict(n=n, lam=lam, max_entry_difference=float(np.max(abs(Hc - M))),
                 nonzeros_referee=int(np.count_nonzero(M)), nonzeros_bond15=int(np.count_nonzero(H)),
                 hermiticity=float(np.max(abs(M - M.conj().T))),
                 charge_commutator=float(np.max(abs(M * (Q[:, None] - Q[None, :])))),
                 )
        out.append(r); print(r)
# Three-site open chain: does the total bad-pair number commute with H?
# The 18-letter parent (arrow grammar) should commute; the 15-letter bond should not.
def three_site(h, pi, dd):
    I = np.eye(dd)
    H = np.kron(h, I) + np.kron(I, h)
    N = np.kron(np.diag(pi), I) + np.kron(I, np.diag(pi))
    return float(np.linalg.norm(H @ N - N @ H, 2))
p = bond18.params(2, K=100., g=.071, td=.003, flux=.34, theta=.23, field=.17, musel=.4)
h18 = bond18.bond(p); _, pen, _, _ = bond18.parts(p); pi18 = (np.diag(pen) > 0).astype(float)
q = my15.default_params(2, K=100., g=.071, td=.003)
pi15 = np.array([float(not my15.allowed(a, b)) for a in my15.LETTERS for b in my15.LETTERS])
r = dict(three_site_commutator_Nbad_parent18=three_site(h18, pi18, 18),
         three_site_commutator_Nbad_15=three_site(my15.bond(q), pi15, 15))
out.append(r); print(r)
json.dump(out, open('check_bond.json', 'w'), indent=1)
