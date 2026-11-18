"""Exhaustive ring checks: every literal transition of the 9- and 10-letter bonds conserves
F = number of (L_a, R_b) adjacencies in the cyclic word with vacancies removed, on ALL d^L words.
So the legal sector F = 0 is exactly invariant at every coupling, with no penalty needed and no vacancy colour."""
import numpy as np, itertools, sys, json, time
import bond9, bond10
def trans(h, d):
    T = {}
    for r, c in zip(*np.nonzero(np.abs(h) > 1e-14)):
        if r != c: T.setdefault(int(c), []).append(int(r))
    return T
def Fcount(word, isL, isR, vac):
    s = [x for x in word if x != vac]
    if len(s) < 2: return 0
    return sum(1 for i in range(len(s)) if isL[s[i]] and isR[s[(i+1)%len(s)]])
def check(mod, d, Lmax, p):
    h = mod.bond(p) if mod is bond10 else (lambda w,o,bad: w)(*bond9.bond(p))
    T = trans(h, d); vac = 9 if d == 10 else -1
    isL = [s in ('L0','L1') for s in mod.LET]; isR = [s in ('R0','R1') for s in mod.LET]
    out = {}
    for L in range(2, Lmax+1):
        viol = ntr = nlegal = 0; t0 = time.time()
        for word in itertools.product(range(d), repeat=L):
            F0 = Fcount(word, isL, isR, vac); nlegal += (F0 == 0)
            for site in range(L):
                j = (site+1) % L
                for r in T.get(word[site]*d + word[j], ()):
                    w2 = list(word); w2[site] = r//d; w2[j] = r % d; ntr += 1
                    if Fcount(w2, isL, isR, vac) != F0: viol += 1
        out[L] = dict(words=d**L, legal_words=nlegal, transitions=ntr, F_violations=viol, seconds=round(time.time()-t0,1))
        print(mod.__name__, 'L=%d'%L, out[L], flush=True)
    return out
if __name__ == '__main__':
    p = dict(lam=0.3,g=0.2,Delta=1.1,D=5.0,nu=0.4,f=0.05,phiA=0.2,td=0.03,J=0.7,t=1.3,mu=0.2,K=3.0)
    res = dict(bond9=check(bond9, 9, int(sys.argv[1]), p), bond10=check(bond10, 10, int(sys.argv[2]), p))
    json.dump(res, open('ring_checks.json','w'), indent=1)
