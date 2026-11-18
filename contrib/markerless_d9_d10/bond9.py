"""Nine-letter markerless route-A bond (folded express rail, no marker letters).
Letters: L0 L1 R0 R1 | W | A B | XA XB.   Cell = L* H R*, H a head.  Only forbidden pair: (L_a, R_b).
All head motion uses one rule, (L_a,H) <-> (H,R_a), so the no-(L,R) sector is exactly invariant.
"""
import numpy as np
LET = ['L0','L1','R0','R1','W','A','B','XA','XB']
IDX = {s:i for i,s in enumerate(LET)}
d = len(LET)
HEADS = ['W','A','B','XA','XB']
NOT_L = ['R0','R1'] + HEADS          # letters that may stand left of a head at the far-left end
NOT_R = ['L0','L1'] + HEADS          # letters that may stand right of a head at the far-right end
CHARGE = {'L0':0,'L1':2,'R0':0,'R1':2,'W':0,'A':-1,'B':1,'XA':-1,'XB':1}

def E(l1,r1,l2,r2):
    m = np.zeros((d*d,d*d),complex); m[IDX[l1]*d+IDX[r1], IDX[l2]*d+IDX[r2]] = 1; return m

def bond(p):
    """returns (w, o, bad): two-site off-diagonal part, on-site energies, projector on forbidden pairs.
    Optional parameters beyond the README's: sw (amplitude of the rail switch, default 1) and phiB (express phase of
    species B, default 0).  The README's bond is sw = 1, phiB = 0."""
    z = {'A':lambda a:(-1)**a, 'B':lambda a:1}
    phi = {'A':p.get('phiA',0.0), 'B':p.get('phiB',0.0)}
    psi = {'A':p.get('psiA',np.pi), 'B':p.get('psiB',np.pi/2)}
    w = np.zeros((d*d,d*d),complex)
    for a in (0,1):
        w += -E(f'L{a}','W','W',f'R{a}')                                  # vacuum walker
        for s in 'AB':
            w += -z[s](a)*E(f'L{a}',s,s,f'R{a}')                          # active token on the data rail
            w += -z[s](a)*np.exp(1j*phi[s])*E(f'L{a}','X'+s,'X'+s,f'R{a}')# blind token on the express rail
    for s in 'AB':
        for y in NOT_L:
            w += -p.get('sw',1.0)*E(y,'X'+s,y,s)                          # rail switch at the far-left end
        for y in NOT_R:
            w += -p.get('td',0.0)*np.exp(1j*psi[s])*E('X'+s,y,s,y)        # weak closing link at the far-right end
    lam = p['lam']
    w += lam*np.exp(-1j*np.pi/4)*E('XA','B','W','W') + lam*np.exp(1j*np.pi/4)*E('XB','A','W','W')   # port
    w += p['g']*E('B','R0','A','R1')                                       # conversion: A + bit1 <-> B + bit0
    w = w + w.conj().T
    o = np.zeros(d)
    nu, f = p.get('nu',0.0), p.get('f',0.0)
    for h in HEADS: o[IDX[h]] = nu
    for s in 'AB':
        o[IDX[s]] += p['Delta']; o[IDX['X'+s]] += p['Delta'] + p['D']
    for s in LET: o[IDX[s]] += f*CHARGE[s]/2
    bad = np.zeros((d*d,d*d))
    for a in (0,1):
        for b in (0,1): bad += E(f'L{a}',f'R{b}',f'L{a}',f'R{b}').real
    return w, o, bad

def full_bond(p):
    """the 81 x 81 bond h9 = w + (1/2)(o x 1 + 1 x o) + K Pi_bad of the README"""
    w, o, bad = bond(p)
    return w + 0.5*(np.kron(np.diag(o), np.eye(d)) + np.kron(np.eye(d), np.diag(o))) + p.get('K', 0.0)*bad

if __name__=='__main__':
    p = dict(lam=0.3,g=0.2,Delta=1.1,D=5.0,nu=0.4,f=0.05,phiA=0.2,td=0.03)
    w,o,bad = bond(p)
    q = np.array([CHARGE[s] for s in LET]); Q2 = np.add.outer(q,q).reshape(-1)
    print('d =',d,' hermitian',np.abs(w-w.conj().T).max(),' charge violation',np.abs(w*(Q2[:,None]-Q2[None,:])).max(),
          ' off-diagonal entries',int((np.abs(w)>0).sum()))
