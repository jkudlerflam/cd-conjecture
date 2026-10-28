"""Literal one-vacancy lift; floating checks use a moderate-scale fixture.

The n=23 infrared hierarchy is NOT tested in floating point by this file.
Run with NumPy and SciPy; output lives beside this script.
"""
from pathlib import Path
import importlib.util
import json
import hashlib
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'contrib/explicit_candidate/bond23.py'
spec = importlib.util.spec_from_file_location('source_bond23', SOURCE)
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def vacancy_lift(h, hopping=1.0, internal_scale=0.001, chemical_potential=0.0):
    """Embed h on occupied pairs and add vacancy swaps; index 0 is vacancy."""
    q = int(round(np.sqrt(h.shape[0])))
    if h.shape != (q*q, q*q):
        raise ValueError('h must have shape (q^2,q^2)')
    d = q+1
    occupied_pairs = np.array([a*d+b for a in range(1,d) for b in range(1,d)])
    out = np.zeros((d*d,d*d),dtype=complex)
    out[np.ix_(occupied_pairs,occupied_pairs)] = internal_scale*h
    for a in range(1,d):
        out[a,a*d] = out[a*d,a] = -hopping
    number = np.r_[0,np.ones(q)]
    out -= chemical_potential/2*np.diag(np.add.outer(number,number).reshape(-1))
    return out


def diagonal_commutator_error(h, values):
    total = np.add.outer(values,values).reshape(-1)
    return float(np.max(np.abs(h*(total[:,None]-total[None,:]))))


def main():
    # Every coupling is moderate here so this tests actual matrix entries.
    p = source.params(td=.013,Delta=2.5,Dx=7.,DM=8.,gp=.031,gm=.027,
                      eps=.017,musel=1.1,K=10.,lam=.023,
                      ph={'0':.21,'+':-.18,'-':.32},
                      phx={'0':.44,'+':.07,'-':-.11},
                      phsw={'0':.12,'+':-.13,'-':.19},
                      phd={'0':.24,'+':.39,'-':-.26})
    h23 = source.bond(p)
    h24 = vacancy_lift(h23,hopping=1.3,internal_scale=.0027,chemical_potential=.0)
    n = np.r_[0,np.ones(23)]
    charge = np.r_[0,[source.charge(a) for a in source.LET]]
    marker = np.r_[0,[int(a in source.MARKERS) for a in source.LET]]
    ix = np.array([a*24+b for a in range(1,24) for b in range(1,24)])
    result = dict(status='numerical_finite_bond_fixture_only',dimension=24,
        bond_shape=list(h24.shape),nonzero_entries=int(np.count_nonzero(h24)),
        hermiticity_error=float(np.max(np.abs(h24-h24.conj().T))),
        occupied_block_error=float(np.max(np.abs(h24[np.ix_(ix,ix)]-.0027*h23))),
        particle_number_error=diagonal_commutator_error(h24,n),
        original_charge_error=diagonal_commutator_error(h24,charge),
        marker_number_error=diagonal_commutator_error(h24,marker),
        alpha_half=.5-2/np.pi**2,
        source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),parameters=p,
        scope='No n=23 CFT, no finite-J charge elimination, no microscopic central charge measured.')
    for k,v in result.items():
        if k.endswith('_error'):
            assert v<1e-13,(k,v)
    np.savez_compressed(HERE/'bond24_fixture.npz',h=h24,n=n,charge=charge,marker=marker)
    (HERE/'bond24_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='parameters'},indent=2))


if __name__=='__main__':
    main()
