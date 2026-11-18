"""Full-sector checks of the complete physical normal-deformation map."""
from __future__ import annotations
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import sys
from pathlib import Path
import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
CONTRIB = HERE.parent
sys.path.insert(0, str(CONTRIB/'qubit_cancellation_classification_20261116'))
import algebra as A
sys.path.insert(0, str(CONTRIB/'qubit_resonant_liquids_20261116'))
import leading

KINETIC = 128/3


def microscopic_series(L, a=0j, zeta=0j, eta=0j, mu=0., beta=0j, pair5=0j, pair6=0j):
    """v=1, g=epsilon chart, alpha=8 conjugate(a), through order six.

    All coefficients are physical onsite/two-site coefficients. Scalars
    merely align the zero energy before each Schur elimination.
    """
    alpha = 8*np.conj(a)
    eye = sparse.eye(1 << L, dtype=complex, format='csr')
    pars = [dict(v=1), dict(a=1, b=-1), dict(u=2, t=-1, p=3),
            dict(p=alpha), dict(u=-34, t=8, p=beta),
            dict(u=-16*np.real(alpha), t=2*np.real(alpha), p=pair5),
            dict(u=196-2*abs(alpha)**2-16*np.real(beta)-KINETIC*mu,
                 t=-80+2*np.real(beta)+KINETIC*zeta, p=pair6, f=KINETIC*eta)]
    hs = [A.chain(L, A.bond(**par)) for par in pars]
    for order, scalar in ((4, 9), (5, 6*np.real(alpha)),
                          (6, abs(alpha)**2+6*np.real(beta)-18)):
        hs[order] = hs[order]+L*scalar*eye
    return hs


def eliminate(L, **parameters):
    h = microscopic_series(L, **parameters)
    states = A.projectors(L)
    first = A.schur_series(h, states)
    keep = np.flatnonzero(abs(first[2].diagonal()) < 1e-9)
    second = A.schur_series(first[2:], keep)
    states = states[keep]
    keep = np.flatnonzero(abs(second[2].diagonal()) < 1e-9)
    third = A.schur_series(second[2:], keep)
    states = states[keep]
    return states, third[2]/KINETIC, dict(
        third_order_norm=A.maxabs(second[1]),
        fifth_order_norm=A.maxabs(third[1]),
        hermiticity_error=A.maxabs(third[2]-third[2].conj().T))


def expected(L, states=None, a=0j, zeta=0j, eta=0j, mu=0., **unused):
    states, h = leading.matrix(L, 3, states=states, alpha=8*np.conj(a),
                               tau=KINETIC*zeta, source=KINETIC*eta)
    numbers = np.array([int(s).bit_count() for s in states])
    return states, h/KINETIC-sparse.diags(mu*numbers)
