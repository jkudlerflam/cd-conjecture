"""Small independent controls for C45's local SW lemma, not the wrapped chain."""
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.linalg import expm, eigvalsh


def kron_all(xs):
    out = np.ones((1, 1))
    for x in xs:
        out = np.kron(out, x)
    return out


def embed(a, support, count, d=4):
    lo, hi = support[0], support[-1] + 1
    assert tuple(support) == tuple(range(lo, hi))
    return np.kron(np.eye(d**lo), np.kron(a, np.eye(d**(count-hi))))


def local_on(a, old, new):
    return embed(a, tuple(i-new[0] for i in old), len(new))


def add_term(terms, support, a):
    terms[support] = terms.get(support, 0) + a


def comm(aa, bb):
    out = {}
    for x, a in aa.items():
        for y, b in bb.items():
            if not set(x).intersection(y):
                continue
            z = tuple(range(min(x[0], y[0]), max(x[-1], y[-1])+1))
            az, bz = local_on(a, x, z), local_on(b, y, z)
            add_term(out, z, az @ bz - bz @ az)
    return out


def plus(*weighted):
    out = {}
    for weight, terms in weighted:
        for x, a in terms.items():
            add_term(out, x, weight*a)
    return out


p1 = np.diag([1., 1., 0., 0.])
a1 = np.diag([0., 0., 1.7, 2.3])
charge = np.array([0, 1, 0, 1])
rng = np.random.default_rng(45)


def conserving_pair():
    q = (charge[:, None] + charge[None, :]).ravel()
    a = rng.normal(size=(16,16)) + 1j*rng.normal(size=(16,16))
    a = (a + a.conj().T)/2
    a[q[:,None] != q[None,:]] = 0
    return a / np.linalg.norm(a, 2)


def split(terms):
    dd, oo, ss = {}, {}, {}
    for x, a in terms.items():
        p = kron_all([p1]*len(x))
        q = np.eye(len(p))-p
        h = sum(embed(a1,(j,),len(x)) for j in range(len(x)))
        energies = h.diagonal()
        inv = np.zeros_like(energies)
        inv[energies > 0] = 1/energies[energies > 0]
        r = np.diag(inv)
        dd[x] = p@a@p + q@a@q
        oo[x] = p@a@q + q@a@p
        ss[x] = r@a@p - p@a@r
    return dd, oo, ss


def full(terms, count):
    return sum(embed(a, x, count) for x, a in terms.items())


vp = conserving_pair()
p2 = np.kron(p1,p1)
vp -= p2@vp@p2
fp = conserving_pair()
rows=[]
for count in (2,3):
    vv = {(j,j+1):vp for j in range(count-1)}
    ff = {(j,j+1):fp for j in range(count-1)}
    vd, vo, s1 = split(vv)
    g2 = plus((1,ff),(1,comm(s1,vd)),(.5,comm(s1,vo)))
    g2d, _, s2 = split(g2)
    a = sum(embed(a1,(j,),count) for j in range(count))
    v, f, t1, t2 = (full(t,count) for t in (vv,ff,s1,s2))
    d1, d2 = full(vd,count), full(g2d,count)
    p = kron_all([p1]*count)
    code = np.flatnonzero(np.diag(p))
    energies = a.diagonal()
    inv = np.zeros_like(energies)
    inv[energies>0] = 1/energies[energies>0]
    expected = p@(f-v@np.diag(inv)@v)@p
    second_error = np.linalg.norm(p@d2@p-expected,2)
    assert second_error < 1e-12
    assert np.linalg.norm(p@d1@p,2) < 1e-12
    ndef = sum(embed(np.eye(4)-p1,(j,),count) for j in range(count))
    vals=[]
    for lam in (.02,.01,.005):
        h = a+lam*v+lam**2*f
        s = lam*t1+lam**2*t2
        u = expm(s)
        d = a+lam*d1+lam**2*d2
        remainder = np.linalg.norm(u@h@u.conj().T-d,2)
        ep = eigvalsh(d[np.ix_(code,code)])[0]
        threshold = eigvalsh(d-ep*np.eye(len(d))-.85*ndef)[0]
        assert threshold > -1e-12
        vals.append({'lambda':lam,'remainder':float(remainder),
                     'remainder_over_lambda3':float(remainder/lam**3),
                     'defect_inequality_min_eigenvalue':float(threshold)})
    rows.append({'cells':count,'dimension':4**count,'code_dimension':2**count,
                 'second_order_identity_error':float(second_error),'runs':vals})

# Exact arithmetic checks of the conservative BCH tail and coupling schedule.
assert Fraction(24,96) < Fraction(1,2)
assert 3*24**2*2*2**3 == 27648
assert 2**16 > 27648 and 2**16 > 144
accepted_checks = 0
for power in range(1,60):
    lam = Fraction(1,2**power)
    # synthetic constants, independent of n; tests algebra, not actual wrapper
    gap, v, G, c, C, eps = map(Fraction,(1,3,20,1,1000000,1))
    eps /= 100
    cap=min(v/G,gap/(192*v),eps*c/(8*C))
    if lam <= cap:
        accepted_checks += 1
        assert C*lam**3/(c*lam**2) <= eps/8
assert accepted_checks > 20
sigma=3*math.sqrt(3/5)-2
assert sigma>0
out={'tag':'numerical: independent local controls, 2 and 3 cells of dimension 4',
     'scope':'Not a wrapper or central-charge calculation.',
     'local_SW':rows,'exact_budget_checks':'passed','irrelevance_margin':sigma}
Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
