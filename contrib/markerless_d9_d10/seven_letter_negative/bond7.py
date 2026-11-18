"""Seven-letter markerless route-A bond.
Letters: L0 L1 R0 R1 W A B.  Cell = L* H R*, H in {W,A,B}.  Only forbidden pair: (L_a, R_b).
"""
import numpy as np, itertools
LET = ['L0','L1','R0','R1','W','A','B']
IDX = {s:i for i,s in enumerate(LET)}
d = 7
def E(l1,r1,l2,r2):
    """|l1 r1><l2 r2| as 49x49"""
    m = np.zeros((d*d,d*d),complex)
    m[IDX[l1]*d+IDX[r1], IDX[l2]*d+IDX[r2]] = 1
    return m
def bond(p):
    """p: dict of parameters. returns (two-site hopping part w, onsite vector o, penalty projector)"""
    w = np.zeros((d*d,d*d),complex)
    ph = {'A':np.exp(1j*p.get('phiA',0)), 'B':np.exp(1j*p.get('phiB',0))}
    z = {'A':lambda a:(-1)**a, 'B':lambda a:1}
    for a in (0,1):
        w += -E(f'L{a}','W','W',f'R{a}')                      # walker: W R_a -> L_a W
        for s in 'AB':
            w += -z[s](a)*ph[s]*E(f'L{a}',s,s,f'R{a}')        # token hop
    lam, chi = p['lam'], p.get('chi',0)
    w += lam*np.exp(1j*chi)*E('A','B','W','W') + lam*np.exp(-1j*chi)*E('B','A','W','W')   # port
    w += p['gR']*E('A','R1','B','R0')                          # right rung: B R0 -> A R1
    w += p.get('gL',0)*E('L1','A','L0','B')                    # left rung:  L0 B -> L1 A
    w = w + w.conj().T
    o = np.zeros(d)
    o[IDX['W']] = p.get('nu',0)
    o[IDX['A']] = p.get('nu',0)+p['DA']; o[IDX['B']] = p.get('nu',0)+p['DB']
    f = p.get('f',0)
    o[IDX['L1']] += f; o[IDX['R1']] += f; o[IDX['A']] += -f/2; o[IDX['B']] += f/2
    bad = np.zeros((d*d,d*d))
    for a in (0,1):
        for b in (0,1):
            bad += E(f'L{a}',f'R{b}',f'L{a}',f'R{b}').real
    return w, o, bad
CHARGE = {'L0':0,'L1':2,'R0':0,'R1':2,'W':0,'A':-1,'B':1}
if __name__=='__main__':
    p = dict(lam=0.3,gR=0.2,gL=0.1*np.exp(0.7j),DA=1.1,DB=1.1,nu=0.4,f=0.05,phiA=0.2)
    w,o,bad = bond(p)
    print('hermitian', np.abs(w-w.conj().T).max())
    q = np.array([CHARGE[s] for s in LET])
    Q2 = np.add.outer(q,q).reshape(-1)
    print('charge violation', np.abs(w*(Q2[:,None]-Q2[None,:])).max())
    print('nonzero offdiag', int((np.abs(w)>0).sum()))
