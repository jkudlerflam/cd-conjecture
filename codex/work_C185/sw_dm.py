"""C185 R4: exact Schrieffer-Wolff (des Cloizeaux) extraction at the Fibonacci point with the transverse DM term
lam (X(x)Z - Z(x)X) added to C173's bond. Built on C179's sw_exact.py conventions (imported unmodified).
Particle = spin down, n = (1 - Z)/2, spin-basis matrix elements.

Predictions (PROOF 3, registered as R4), with Delta = 0, w = 0:
  NN hop  <j|H|j+1>  = +4 lam^2/V
  on-site E(j)-E(0)+mu = -8 lam^2/V
  V2 = +8 lam^2/V,  V3 = 0,  NNN hop = 0
With Delta != 0: the Delta terms of C173 add, and dN = 1 elements of order Delta lam / V appear.
"""
import numpy as np, sys, importlib.util
spec = importlib.util.spec_from_file_location("sw", "../work_C179/sw_exact.py")
sw = importlib.util.module_from_spec(spec); spec.loader.exec_module(sw)
op, X, Y, Z = sw.op, sw.X, sw.Y, sw.Z

def heff_dm(L, Delta, V, w, mu, lam, lamy=0.0):
    Jz = V/4; hz = (mu - 4*Jz)/2
    Jx = (w.real + Delta)/2; Jy = (w.real - Delta)/2; D = w.imag/2
    H = sw.full_chain(L, Jx, Jy, D, Jz, hz)
    H += (hz/2)*(op(L, {0: Z}) + op(L, {L-1: Z})) + Jz*(op(L, {0: Z}) + op(L, {L-1: Z}))
    for j in range(L-1):
        H += lam*(op(L, {j: X, j+1: Z}) - op(L, {j: Z, j+1: X}))
        H += lamy*(op(L, {j: Y, j+1: Z}) - op(L, {j: Z, j+1: Y}))
    e, v = np.linalg.eigh(H)
    occ = lambda s: [(s >> (L-1-j)) & 1 for j in range(L)]
    Pidx = [s for s in range(2**L) if not any(a & b for a, b in zip(occ(s)[:-1], occ(s)[1:]))]
    m = len(Pidx)
    A = v[Pidx, :m]
    u, sv, vh = np.linalg.svd(A)
    He = (u @ vh) @ np.diag(e[:m]) @ (u @ vh).conj().T
    return He, {s: i for i, s in enumerate(Pidx)}, e[m]-e[m-1], sv.min()

if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    c = L//2 - 1
    for Delta in (0.0, 1.0):
        print(f"\n=== Delta = {Delta}, w = 0, L = {L}. lam = 1 when Delta = 0, lam = 0.25 Delta otherwise ===")
        lam = 1.0 if Delta == 0 else 0.25
        for V in (40.0, 80.0, 160.0, 320.0):
            mu = 0.0
            He, pos, gap, smin = heff_dm(L, Delta, V, 0j, mu, lam)
            el = lambda f, i: sw.elem(He, pos, L, f, i)
            E0 = el([], []); Ej = lambda *s: el(list(s), list(s))
            u = lam**2/V; t2 = Delta**2/(2*V)
            rows = [("hop1", el([c], [c+1]), 4*u),
                    ("hop2", el([c], [c+2]), -t2),
                    ("onsite", (Ej(c)-E0).real, -8*u + 3*Delta**2/V),
                    ("V2", (Ej(c, c+2)-Ej(c)-Ej(c+2)+E0).real, 8*u - Delta**2/V),
                    ("V3", (Ej(c-1, c+2)-Ej(c-1)-Ej(c+2)+E0).real, -Delta**2/(3*V)),
                    ("V4", (Ej(c-1, c+3)-Ej(c-1)-Ej(c+3)+E0).real, 0.0),
                    ("dN1 <c|H|0>", el([c], []), None),
                    ("dN1 <c,c+2|H|c+2>", el([c, c+2], [c+2]), None),
                    ("dN1 <c,c+3|H|c+3>", el([c, c+3], [c+3]), None),
                    ("dN1 <c,c+4|H|c+4>", el([c, c+4], [c+4]), None)]
            print(f" V = {V}: gap {gap:.2f}, smin {smin:.5f}")
            for nm, me, pr in rows:
                if pr is None:
                    print(f"   {nm:20s} {complex(me):.5e}   x V/(Delta lam) = {complex(me)*V/max(Delta*lam,1e-30):.4f}")
                else:
                    r = complex(me)/pr if abs(pr) > 1e-15 else float('nan')
                    print(f"   {nm:20s} {complex(me):.5e}  pred {pr:.5e}  ratio {complex(r):.4f}")
