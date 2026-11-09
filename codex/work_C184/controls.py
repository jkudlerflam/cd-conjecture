"""Controls for the two-size estimator at the same sizes: XX (c=1), critical TFI (c=1/2) as qubit bonds, and on d = 4 rungs
two decoupled XX chains (c = 2) and TFI x XX (c = 3/2). Usage: python3 controls.py L chis out.jsonl label"""
import sys, json, time
import numpy as np
from dmrg_d import dmrg
X = np.array([[0,1],[1,0]],float); Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1.,-1.]); I = np.eye(2)
k = np.kron
def xx(): return (k(X,X) + k(Y,Y)).real
def tfi(): return -k(X,X) - 0.5*(k(Z,I) + k(I,Z))
def rung(hA, hB):  # sites (a1 b1)(a2 b2), d=4; hA acts on a1 a2, hB on b1 b2
    A = hA.reshape(2,2,2,2); B = hB.reshape(2,2,2,2); Ia = np.eye(2)
    H = np.einsum('ikjl,mn,op->imkojnlp', A, Ia, Ia).reshape(16,16) * 0
    # build explicitly: basis index (a1,b1,a2,b2)
    H = np.zeros((2,2,2,2,2,2,2,2))
    H += np.einsum('acAC,bB,dD->abcdABCD', hA.reshape(2,2,2,2), Ia, Ia)
    H += np.einsum('bdBD,aA,cC->abcdABCD', hB.reshape(2,2,2,2), Ia, Ia)
    return H.reshape(16,16)
models = {'xx': xx(), 'tfi': tfi(), 'xx2': rung(xx(), xx()), 'tfixx': rung(tfi(), xx())}
if __name__ == '__main__':
    L = int(sys.argv[1]); chis = tuple(int(c) for c in sys.argv[2].split(',')); fout = sys.argv[3]; lab = sys.argv[4]
    h = models[lab]; t0 = time.time()
    if lab == 'tfi':   # restore full field on edge sites
        ons = [None]*L; ons[0] = -0.5*Z; ons[L-1] = -0.5*Z
        r = dmrg(h, L, chis=chis, onsite=ons)
    else:
        r = dmrg(h, L, chis=chis)
    rec = dict(label='ctrl_'+lab, L=L, chis=list(chis), E=r['E'], trunc=r['hist'][-1][2], dE_last=r['hist'][-1][1]-r['hist'][-2][1],
               S_mid=float(r['S'][L//2-1]), S=[round(float(s),6) for s in r['S']], sec=round(time.time()-t0,1))
    open(fout,'a').write(json.dumps(rec)+'\n'); print(lab, L, rec['S_mid'], rec['trunc'], rec['sec'])
