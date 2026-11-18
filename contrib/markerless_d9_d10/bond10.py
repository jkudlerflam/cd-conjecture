"""Ten-letter mobile bond: the nine material letters of bond9 plus ONE uncoloured, fully transparent vacancy V.

    h10 = J * embed(b9) - t * sum_a (|a V><V a| + h.c.) - (mu/2)(n_occ x 1 + 1 x n_occ),
    b9  = w + (1/2)(o x 1 + 1 x o) + K * P(L_a R_b adjacent)        (the frozen bond without its field convention changed)

embed() acts on material-material pairs only, so the material on-site energies are contact terms, exactly as in
C88 (contrib/small_d_note, eq. h13).  Then a fast vacancy fluid multiplies the WHOLE material bond by the contact
probability and the frozen selector and Bloch data are unchanged to first order in J/t.

`onsite='genuine'` gives the first version of this bond, in which the on-site energies did not need contact.
The referee showed that it does not select the cell pattern once vacancies are present (REFEREE.md, section 5);
it is kept only so that mobile_selector.py can reproduce that failure.
"""
import numpy as np
import bond9
LET = bond9.LET + ['V']; IDX = {s:i for i,s in enumerate(LET)}; d = 10
CHARGE = dict(bond9.CHARGE, V=0)
def bond(p, onsite='contact'):
    w9, o9, bad9 = bond9.bond(p)
    J, t, mu, K = p.get('J',1.0), p.get('t',1.0), p.get('mu',0.0), p.get('K',0.0)
    m = [i*d+j for i in range(9) for j in range(9)]            # material-material pair indices in the 10-letter basis
    emb = np.zeros((d*d, d*d), complex)
    b9 = w9 + K*bad9
    if onsite == 'contact':
        b9 = b9 + 0.5*(np.kron(np.diag(o9), np.eye(9)) + np.kron(np.eye(9), np.diag(o9)))
    emb[np.ix_(m, m)] = J*b9
    o = np.zeros(d); o[:9] = -mu
    if onsite == 'genuine': o[:9] += J*o9
    onsite_term = 0.5*(np.kron(np.diag(o), np.eye(d)) + np.kron(np.eye(d), np.diag(o)))
    hop = np.zeros((d*d, d*d))
    for a in range(9):
        hop[a*d+9, 9*d+a] += -t; hop[9*d+a, a*d+9] += -t          # |a V><V a| + h.c.
    return emb + onsite_term + hop
if __name__ == '__main__':
    p = dict(lam=0.3,g=0.2,Delta=1.1,D=5.0,nu=0.4,f=0.05,phiA=0.2,td=0.03,J=0.7,t=1.3,mu=0.2,K=3.0)
    for conv in ('contact', 'genuine'):
        h = bond(p, conv)
        q = np.array([CHARGE[s] for s in LET]); Q2 = np.add.outer(q,q).reshape(-1)
        occ = np.array([1]*9+[0]); N2 = np.add.outer(occ,occ).reshape(-1)
        hd = np.array([1 if s in bond9.HEADS else 0 for s in LET]); H2 = np.add.outer(hd,hd).reshape(-1)
        print(conv, 'd=10 hermitian', np.abs(h-h.conj().T).max(), ' Q', np.abs(h*(Q2[:,None]-Q2[None,:])).max(),
              ' N_occ', np.abs(h*(N2[:,None]-N2[None,:])).max(), ' N_head', np.abs(h*(H2[:,None]-H2[None,:])).max(),
              ' offdiag', int((np.abs(h-np.diag(np.diag(h)))>1e-14).sum()))
