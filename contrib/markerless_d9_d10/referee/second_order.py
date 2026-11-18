"""Item 3: second-order (in the port) logical Hamiltonian from the literal letter model.

Method (referee's own):  H = H0 + g C + lam V  is built literally from the referee's
h9 on a chain of sites.  P = code space (clock ground state of every cell, any bits).
Because C annihilates P and one conversion leaves two tokens of equal species,

    M0 = - P V R V P                (lam^2 g^0)
    M1 = + P V R C R V P            (lam^2 g^1)
    M2 = - P V R C R C R V P        (lam^2 g^2, exact coefficient)

with R = (H0 - E0)^{-1} on the token-pair sector, computed by sparse LU.
No resolvent factorisation over cells is assumed.

Systems:
   'pad'   : [padL] cell cell [padR]   open chain (author's geometry is padL=R0, padR=L0)
   'ring3' : three cells on a periodic ring (no padding)
   'ring2' : two cells on a periodic ring (shows the wrap-around artefact)
"""
import sys, os, json, itertools, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_bonds as rb
import chain as ch

L_OF = (rb.L0, rb.L1)
R_OF = (rb.R0, rb.R1)


def cell_word(bits, j, head=rb.W):
    return tuple(L_OF[a] for a in bits[:j]) + (head,) + tuple(R_OF[a] for a in bits[j:])


def clock_vec(n):
    j = np.arange(n + 1)
    return np.sqrt(2.0 / (n + 2)) * np.sin(np.pi * (j + 1) / (n + 2))


def par_of(n, a, D, td, lam, g, nu=0.0, f=0.0, K=0.0, phiA=None, psiA=np.pi, psiB=np.pi / 2):
    return dict(nu=nu, Delta=a - 2 * np.cos(np.pi / (n + 2)), D=D, f=f, K=K, td=td, g=g, lam=lam,
                phiA=np.pi / (2 * (n + 1)) if phiA is None else phiA, psiA=psiA, psiB=psiB)


def build(n, ncell, geometry, base, padL=rb.R0, padR=rb.L0):
    """returns dict with H0, C, V (sparse), the code isometry Psi (columns = bit strings), E0"""
    periodic = geometry.startswith("ring")
    bitsets = list(itertools.product((0, 1), repeat=n))
    v = clock_vec(n)
    seeds = []
    code_entries = []   # (column, word, amplitude)
    cols = list(itertools.product(range(len(bitsets)), repeat=ncell))
    for c, combo in enumerate(cols):
        for js in itertools.product(range(n + 1), repeat=ncell):
            wd = ()
            amp = 1.0
            for cell in range(ncell):
                wd += cell_word(bitsets[combo[cell]], js[cell])
                amp *= v[js[cell]]
            if not periodic:
                wd = (padL,) + wd + (padR,)
            seeds.append(wd)
            code_entries.append((c, wd, amp))
    generic = dict(base); generic.update(lam=0.37, g=0.23, td=base["td"] if base["td"] != 0 else 0.0)
    hg = rb.h9(generic)
    words, index = ch.closure(seeds, hg, 9, periodic)
    def Hof(lam, g):
        p = dict(base); p.update(lam=lam, g=g)
        H, leak = ch.hamiltonian(words, index, rb.h9(p), 9, periodic)
        assert leak == 0
        return H
    H0 = Hof(0.0, 0.0)
    C = Hof(0.0, 1.0) - H0
    V = Hof(1.0, 0.0) - H0
    # check linearity in (lam, g): H(lam,g) = H0 + lam V + g C
    lin = abs(Hof(0.3, 0.7) - (H0 + 0.3 * V + 0.7 * C)).max()
    rows = [index[wd] for (_, wd, _) in code_entries]
    Psi = sp.coo_matrix(([e[2] for e in code_entries], (rows, [e[0] for e in code_entries])),
                        shape=(len(words), len(cols)), dtype=complex).tocsc()
    HPsi = H0 @ Psi
    E0 = (Psi.conj().T @ HPsi).diagonal().real
    resid = abs(HPsi - Psi @ sp.diags(E0)).max()
    # sectors
    letters = np.array(words)
    ntok = ((letters == rb.A) | (letters == rb.B) | (letters == rb.XA) | (letters == rb.XB)).sum(1)
    nhead = ntok + (letters == rb.W).sum(1)
    nLR = np.zeros(len(words), int)
    N = letters.shape[1]
    for i in range(N - 1 + (1 if periodic else 0)):
        j = (i + 1) % N
        nLR += (np.isin(letters[:, i], L_OF) & np.isin(letters[:, j], R_OF))
    return dict(n=n, ncell=ncell, words=words, index=index, H0=H0, C=C, V=V, Psi=Psi, E0=E0,
                E0_spread=float(E0.max() - E0.min()), code_resid=float(resid), linearity=float(lin),
                ntok=ntok, nhead=nhead, F_max=int(nLR.max()), dim=len(words), bitsets=bitsets, cols=cols)


def effective(sys_):
    H0, C, V, Psi = sys_["H0"], sys_["C"], sys_["V"], sys_["Psi"]
    E0 = float(sys_["E0"].mean())
    ntok = sys_["ntok"]
    tok = np.nonzero(ntok > 0)[0]
    vac = np.nonzero(ntok == 0)[0]
    VP = (V @ Psi)
    first = (Psi.conj().T @ VP).toarray()
    # V P must live entirely in the token sector
    vp_vac = abs(VP[vac, :]).max() if len(vac) else 0.0
    A = (H0[tok, :][:, tok] - E0 * sp.identity(len(tok), dtype=complex, format="csc")).tocsc()
    lu = spla.splu(A)
    Ct = C[tok, :][:, tok].tocsc()
    b = VP[tok, :].toarray()
    X1 = lu.solve(b)
    X2 = lu.solve(Ct @ X1)
    X3 = lu.solve(Ct @ X2)
    M0 = -(b.conj().T @ X1)
    M1 = +(b.conj().T @ X2)
    M2 = -(X1.conj().T @ (Ct @ X2))
    M4 = -(X2.conj().T @ (Ct @ X3))     # lam^2 g^4 coefficient, for the size of the tail
    # smallest eigenvalue of the token sector above E0 (denominator scale)
    return dict(first=first, vp_vac=float(vp_vac), M0=M0, M1=M1, M2=M2, M4=M4, E0=E0)


def fermion_ops(nmodes):
    """Jordan-Wigner c_m on (C^2)^{nmodes}; basis index = bits as big-endian integer
    (mode 0 most significant), matching itertools.product order."""
    I = np.eye(2); Z = np.diag([1.0, -1.0]); sm = np.array([[0, 1], [0, 0]], float)  # |0><1|
    ops = []
    for m in range(nmodes):
        mats = [Z] * m + [sm] + [I] * (nmodes - m - 1)
        out = np.array([[1.0]])
        for x in mats:
            out = np.kron(out, x)
        ops.append(out)
    return ops


def one_body(M, nmodes):
    """extract const and h from the 0- and 1-particle sectors, and the Gaussian residual"""
    dim = 2 ** nmodes
    c = fermion_ops(nmodes)
    vac = np.zeros(dim); vac[0] = 1.0
    const = vac @ M @ vac
    states = [c[m].T @ vac for m in range(nmodes)]       # c_m^dag |0>
    h = np.array([[states[m].conj() @ M @ states[mp_] for mp_ in range(nmodes)] for m in range(nmodes)]) \
        - const * np.eye(nmodes)
    G = const * np.eye(dim, dtype=complex)
    for m in range(nmodes):
        for mp_ in range(nmodes):
            if h[m, mp_] != 0:
                G += h[m, mp_] * (c[m].T @ c[mp_])
    return const, h, float(abs(M - G).max())


def analyse(n, a, D, td, geometry="pad", ncell=2, padL=rb.R0, padR=rb.L0, nu=0.3, f=0.0, K=0.0, **kw):
    base = par_of(n, a, D, td, 0.0, 0.0, nu=nu, f=f, K=K, **kw)
    s = build(n, ncell, geometry, base, padL, padR)
    e = effective(s)
    nm = n * ncell
    const0 = np.trace(e["M0"]) / e["M0"].shape[0]
    c2, h, gauss = one_body(e["M2"], nm)
    out = dict(n=n, a=a, D=D, td=td, geometry=geometry, ncell=ncell, dim=s["dim"], F_max=s["F_max"],
               code_resid=s["code_resid"], E0_spread=s["E0_spread"], linearity=s["linearity"],
               first_order=float(abs(e["first"]).max()), VP_in_vacuum_sector=e["vp_vac"],
               M0_const=complex(const0), M0_bit_dependence=float(abs(e["M0"] - const0 * np.eye(len(e["M0"]))).max()),
               M1_max=float(abs(e["M1"]).max()), M2_max=float(abs(e["M2"]).max()),
               M2_hermiticity=float(abs(e["M2"] - e["M2"].conj().T).max()),
               M2_const=complex(c2), gauss_resid=gauss, hop_scale=float(abs(h).max()),
               M4_max=float(abs(e["M4"]).max()))
    return out, h, e, s


if __name__ == "__main__":
    np.set_printoptions(precision=6, linewidth=200, suppress=False)
    out, h, e, s = analyse(2, 2.6, 3.0, 0.07)
    print(json.dumps({k: (str(v) if isinstance(v, complex) else v) for k, v in out.items()}, indent=1))
    print("one-body matrix h (units lam^2 g^2):\n", h)
