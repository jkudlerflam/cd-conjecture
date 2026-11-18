"""Exact two-body/bare-band checks of a minimal physical real-hop detuning.

No Gaussian parameters are assigned to these band numbers. The pair bound
state is a vacuum obstruction, not a finite-density CFT calculation.
"""
from quartet import np, la, threadpool_limits
from scipy.optimize import brentq
from pathlib import Path
import json
import time


def bound_pair(kappa):
    # Product q=z1*z2 of the two decaying relative-coordinate roots.
    q=brentq(lambda q:(1-q)*(3-q)*(1+q)**2-kappa**2*q**3,1e-12,1,
             xtol=1e-14)
    p=-kappa*q/(1+q)
    E=-2*(q+1/q+kappa**2*q/(1+q)**2)
    edge=-4-kappa**2/2
    roots=np.roots([1,-p,q])
    assert max(abs(roots))<1
    # Independent literal relative-coordinate half-line calculation.
    size=512
    band=np.zeros((3,size)); band[0,0]=-2; band[0,1]=-2/3
    band[1,:-1]=2*kappa; band[2,:-2]=2
    box=float(la.eig_banded(band,lower=True,eigvals_only=True,
                            select='i',select_range=(0,0))[0])
    assert abs(E-box)<2e-8
    return dict(kappa=kappa,q=q,bound_energy=E,continuum_edge=edge,
                binding_energy=edge-E,relative_box=512,
                box_error=abs(E-box),decay_modulus=float(max(abs(roots))))


def band_minima(kappa,real_hop):
    derivative=lambda k:-2*real_hop*np.sin(k)+2*kappa*np.cos(k)+4*np.sin(2*k)
    dispersion=lambda k:2*real_hop*np.cos(k)+2*kappa*np.sin(k)-2*np.cos(2*k)
    # Two known isolated wells; bounded scalar roots, no parameter scan.
    ks=[brentq(derivative,-np.pi-.8,-np.pi+.8),
        brentq(derivative,-.8,.8)]
    return np.array(ks),np.array([dispersion(k) for k in ks])


def calculate():
    bindings=[bound_pair(k) for k in (.5,1,2)]
    # Physical coefficients of the five-monomial detuned example. The field
    # and every coefficient are fixed independently of microscopic length.
    V=16; a=33/64; b=31/64; u=4; h=127/16; v=1/512
    Delta=np.sqrt((a+b)**2+v*v)
    t=Delta**2/(2*V)
    kappa=(a-b)/t; real=v/t
    mu=(V-2*h-3*Delta**2/V)/t
    ks,es=band_minima(kappa,real)
    q4=float(2*(sum(ks)+np.pi))
    baseks,basees=band_minima(1,1/16)
    # At total momentum pi the real nearest-neighbour hop cancels exactly.
    pair=bound_pair(kappa)
    return dict(bindings=bindings,
                detuned_example=dict(coefficients=dict(XY=a,YX=b,ZZ=u,Z=-h,XX=v),
                    physical_dimension=2,range=1,one_site_translation=True,
                    independent_real_coefficients_before_scale=5,
                    ratios_after_scale=4,literal_Pauli_monomials=5,
                    literal_nonzero_spread=h/v,
                    Delta=Delta,t=t,kappa=kappa,real_hop_over_t=real,
                    chemical_potential_over_t=mu,
                    bare_valley_centres=ks.tolist(),bare_minimum_energies_over_t=es.tolist(),
                    bare_quartet_mismatch=q4,bare_valley_offset=float(abs(es[0]-es[1])),
                    both_bare_wells_below_mu=bool(max(es)<mu),
                    pi_bound_pair_above_two_lowest_atoms=float(pair['bound_energy']-2*min(es)),
                    leading_separation2_source_over_t=-2*(Delta/V)*real,
                    leading_separation3_source_over_t=4*(Delta/V)/3,
                    scope='Microscopic detuning with tested leading kinematics; no stable phase established.'),
                linearized_at_kappa1=dict(real_hop=1/16,
                    predicted_quartet_mismatch=-4*(1/16)/15,
                    actual_quartet_mismatch=float(2*(sum(baseks)+np.pi)),
                    predicted_valley_offset=float(np.sqrt(15)/16),
                    actual_valley_offset=float(abs(basees[0]-basees[1]))))


if __name__=='__main__':
    start=time.perf_counter()
    with threadpool_limits(limits=1): out=calculate()
    out['seconds']=time.perf_counter()-start
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
