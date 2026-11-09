"""C186: C175's dilute pseudospin criterion (C175 PROOF 1, scat.py imported read-only from work_C175) for the
number-conserving cusp model with V2 = -2 + v, V3 = -2/3, past the pair boundary. c = 2 iff a_- < a_s < a_+,
with a_+- the eigenvalues of the Hermitian part of the unlike threshold matrix A and a_s the like length."""
import sys, numpy as np
sys.path.insert(0, "../work_C175")
from scat import A_matrix, a_like
for u in (4*np.sin(np.pi/12), 1.0, 2.0):
    vc = u*u/4
    for f in (0.0, 0.8, 1.2, 1.5, 2.0, 3.0):
        v = f*vc
        A, ka, kb, K = A_matrix(1j*u, V2=-2.0+v, V3=-2/3)
        AH = (A + A.conj().T)/2
        am, ap = np.linalg.eigvalsh(AH)
        a_s, _ = a_like(1j*u, 0, V2=-2.0+v, V3=-2/3)
        Ad, Ao = AH[0, 0].real, abs(AH[0, 1])
        ds = (a_s.real - Ad)/Ao
        verdict = "c=2 window" if am < a_s.real < ap else ("vector chiral" if a_s.real <= am else "valley crystal")
        print(f"u={u:.4f} v={v:.4f} (v/v_c={f}) a_-={am:+.3f} a_+={ap:+.3f} a_s={a_s.real:+.3f} Delta_sigma={ds:+.3f} -> {verdict}"
              + ("   [pair bound: criterion not applicable]" if f < 1 else ""), flush=True)
