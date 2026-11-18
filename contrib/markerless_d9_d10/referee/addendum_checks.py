"""Addendum: checks of the corrected mobile bond (contact on-site terms), referee's own code.

  h10 = J embed(b9) - t sum_a (|a V><V a| + h.c.) - (mu/2)(n_occ x 1 + 1 x n_occ),
  b9  = w + (1/2)(o x 1 + 1 x o) + K Pi_bad,  embedded on material-material pairs only.

 A. entry-by-entry comparison with bond10.bond(p) at random generic points
 B. ring invariants (F, Q, heads, occupied, squeezed boundary pattern) for L <= 5
 C. one cell + Nv vacancies against the headless ring: eps_cell, its ratio to the frozen value,
    and the contact probability c of the vacancy fluid measured in the headless ground state
 D. bit dependence at zero port with a nonzero field f (the field is now a contact term too)
"""
import sys, os, json, itertools
import numpy as np
import scipy.sparse.linalg as spla
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.dirname(HERE))
import ref_bonds as rb, chain as ch, ring_invariants as ri

def h10_contact(par):
    b9 = rb.h9(par)                                   # w + (1/2)(o x 1 + 1 x o) + K Pi_bad
    H = np.zeros((100, 100), dtype=complex)
    for i in range(9):
        for j in range(9):
            for k in range(9):
                for l in range(9):
                    H[i * 10 + j, k * 10 + l] = par["J"] * b9[i * 9 + j, k * 9 + l]
    nocc = np.diag([1.0] * 9 + [0.0]); I = np.eye(10)
    H += -0.5 * par["mu"] * (np.kron(nocc, I) + np.kron(I, nocc))
    for a in range(9):
        H[a * 10 + rb.V, rb.V * 10 + a] += -par["t"]; H[rb.V * 10 + a, a * 10 + rb.V] += -par["t"]
    return H

out = {}
# ---- A
import bond10
rng = np.random.default_rng(20261113)
dmax = 0.0; dgen = 0.0
for trial in range(12):
    par = rb.random_par(rng, mobile=True)
    pa = dict(lam=par["lam"], g=par["g"], Delta=par["Delta"], D=par["D"], nu=par["nu"], f=par["f"], phiA=par["phiA"],
              psiA=par["psiA"], psiB=par["psiB"], td=par["td"], J=par["J"], t=par["t"], mu=par["mu"], K=par["K"])
    dmax = max(dmax, float(abs(bond10.bond(pa) - h10_contact(par)).max()))
    old = dict(par); old["K"] = par["K"]      # first version, with the README's J*K now also in the code?
    dgen = max(dgen, float(abs(bond10.bond(pa, onsite="genuine") - rb.h10(old)).max()))
Hc = h10_contact(par)
out["A"] = dict(points=12, max_entry_difference_contact=dmax, max_entry_difference_genuine_vs_referee_first_version=dgen,
                hermiticity=float(abs(Hc - Hc.conj().T).max()),
                offdiag_entries=int((abs(Hc - np.diag(np.diag(Hc))) > 1e-14).sum()))
print(json.dumps(out["A"]), flush=True)
# ---- B
out["B"] = []
for L in range(2, 6):
    r = ri.run(10, L, Hc, ["F", "Q", "heads", "occ", "sqz"])
    out["B"].append(r); print(json.dumps(r), flush=True)

# ---- C
def nu_sel(b):
    a = lambda x: 2 * (1 - np.cos(np.pi / (x + 1)))
    dl = lambda x: (x + 1) * a(x) - x * a(x + 1)
    return 2 - (dl(b) + dl(b - 1)) / 2

def placements(material, Nv):
    M = len(material); L = M + Nv
    res = []
    for vac in itertools.combinations(range(L), Nv):
        w = [rb.V] * L; it = iter(material)
        for i in range(L):
            if i not in vac: w[i] = next(it)
        res.append(tuple(w))
    return res

def ground(seeds, par, want_vec=False):
    words, index = ch.closure(seeds, h10_contact(par), 10, True)
    H, leak = ch.hamiltonian(words, index, h10_contact(par), 10, True)
    assert leak == 0
    if H.shape[0] <= 600:
        w, v = np.linalg.eigh(H.toarray()); E, vec = float(w[0]), v[:, 0]
    else:
        w, v = spla.eigsh(H, k=1, which="SA", tol=1e-12); E, vec = float(w[0]), v[:, 0]
    return (E, words, vec) if want_vec else E

def contact_probability(words, vec, M):
    L = len(words[0]); tot = 0.0
    for wd, amp in zip(words, vec):
        nn = sum(1 for i in range(L) if wd[i] != rb.V and wd[(i + 1) % L] != rb.V)
        tot += abs(amp) ** 2 * nn
    return tot / M

out["C"] = []
for n, Nvmax in ((4, 6), (9, 4)):
    b = n + 1; nu = nu_sel(b); frozen = nu - 2 * np.cos(np.pi / (b + 1))
    for t in (25.0, 5.0, 1.0, 0.2):
        par = dict(nu=nu, Delta=0.5, D=3.0, f=0.0, K=2.0, td=0.0, g=0.0, lam=0.0, phiA=0.1, psiA=np.pi, psiB=np.pi / 2, J=1.0, t=t, mu=0.0)
        for Nv in range(0, Nvmax + 1):
            EA = ground(placements((rb.W,) + (rb.R0,) * n, Nv), par)
            EB, wB, vB = ground(placements((rb.R0,) * b, Nv), par, want_vec=True)
            c = contact_probability(wB, vB, b)
            d = dict(n=n, t_over_J=t, Nv=Nv, density=Nv / (b + Nv), eps_cell=EA - EB, frozen=frozen, ratio=(EA - EB) / frozen, contact_probability=c)
            out["C"].append(d); print(json.dumps(d), flush=True)

# ---- D: bit dependence at zero port with a field, one n = 3 cell with a single 1-bit at position j
out["D"] = []
n, Nv = 3, 3; b = n + 1; nu = nu_sel(b)
for t in (25.0, 1.0):
    for f in (0.0, 1e-3):
        par = dict(nu=nu, Delta=0.5, D=3.0, f=f, K=2.0, td=0.0, g=0.0, lam=0.0, phiA=0.1, psiA=np.pi, psiB=np.pi / 2, J=1.0, t=t, mu=0.0)
        E = {}
        for bits in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]:
            cell = (rb.W,) + tuple(rb.R1 if a else rb.R0 for a in bits)
            E[bits] = ground(placements(cell, Nv), par)
        _, wB, vB = ground(placements((rb.R0,) * b, Nv), par, want_vec=True)
        c = contact_probability(wB, vB, b)
        shifts = [E[bt] - E[(0, 0, 0)] for bt in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
        d = dict(n=n, Nv=Nv, t_over_J=t, f=f, contact_probability=c, bit_energy_shift=shifts,
                 shift_over_f=[s / f for s in shifts] if f else None,
                 spread_over_f=(max(shifts) - min(shifts)) / f if f else None)
        out["D"].append(d); print(json.dumps(d), flush=True)
json.dump(out, open(os.path.join(HERE, "addendum_checks.json"), "w"), indent=1)
