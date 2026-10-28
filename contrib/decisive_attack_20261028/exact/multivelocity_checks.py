"""Exact-formula numerical checks for MULTIVELOCITY_REVIEW.md; numpy only."""
import json
from pathlib import Path
import numpy as np

rng = np.random.default_rng(2026102807)
records = []
for branches in (2, 4, 6, 8):
    velocity = np.linspace(.4, 1.7, branches)*np.where(np.arange(branches)%2, -1, 1)
    v = np.diag(velocity)
    sign = np.diag(np.sign(velocity))
    d = np.diag(1/(4*np.pi*np.abs(velocity)))
    sqrt_d = np.diag(np.sqrt(np.diag(d)))
    raw = rng.normal(size=(branches, branches))
    coupling = .3*(raw+raw.T)/2
    enorm = np.linalg.norm(sqrt_d @ coupling @ sqrt_d, 2)
    coupling *= .25/max(enorm, .25)
    enorm = np.linalg.norm(sqrt_d @ coupling @ sqrt_d, 2)
    identity = np.eye(branches)
    m = identity-d@coupling
    nmat = identity+d@coupling
    minv = np.linalg.inv(m)
    k = minv @ sign @ minv.T/(2*np.pi)
    h = (np.linalg.inv(d)-coupling@d@coupling)/2
    a = minv @ v @ nmat
    kh_error = float(np.max(np.abs(a-k@h)))
    he, ho = np.linalg.eigh(h)
    assert he.min() > 0 and enorm < 1
    sqrt_h = (ho*np.sqrt(he)) @ ho.T
    speeds, modes = np.linalg.eigh(sqrt_h @ k @ sqrt_h)
    assert np.sum(speeds > 0) == np.sum(velocity > 0)
    assert np.sum(speeds < 0) == np.sum(velocity < 0)
    transform = np.diag(1/np.sqrt(2*np.pi*np.abs(speeds))) @ modes.T @ sqrt_h
    z = np.diag(1+.1*rng.random(branches))
    source = identity+.03*rng.normal(size=(branches, branches))
    physical_mix = source @ np.linalg.inv(z)
    physical_transform = transform @ np.linalg.inv(physical_mix)
    max_contact_error = max_mode_error = max_physical_error = 0.0
    for omega, momentum in ((.9, .3), (.6, -.8), (1.3, 0), (.01, 1.2)):
        bubble = d @ np.diag((-1j*omega+velocity*momentum)/(1j*omega+velocity*momentum))
        covariance = np.linalg.solve(identity+bubble@coupling, bubble)
        contact = -minv@d
        subtracted = covariance-contact
        rational = momentum*np.linalg.solve(1j*omega*identity+a*momentum, k)
        mode_covariance = transform @ subtracted @ transform.T
        target = np.diag(np.sign(speeds)*momentum/(2*np.pi*(1j*omega+speeds*momentum)))
        physical_covariance = physical_mix @ subtracted @ physical_mix.T
        physical_mode = physical_transform @ physical_covariance @ physical_transform.T
        max_contact_error = max(max_contact_error, float(np.max(np.abs(subtracted-rational))))
        max_mode_error = max(max_mode_error, float(np.max(np.abs(mode_covariance-target))))
        max_physical_error = max(max_physical_error, float(np.max(np.abs(physical_mode-target))))
    assert max(kh_error, max_contact_error, max_mode_error, max_physical_error) < 1e-11
    records.append(dict(branches=branches, stability_norm=float(enorm),
                        smallest_energy_eigenvalue=float(he.min()),
                        smallest_absolute_speed=float(np.min(np.abs(speeds))),
                        positive_speeds=int(np.sum(speeds > 0)),
                        negative_speeds=int(np.sum(speeds < 0)),
                        kh_identity_error=kh_error,
                        contact_identity_error=max_contact_error,
                        diagonal_mode_error=max_mode_error,
                        physical_source_mode_error=max_physical_error))
# Equal velocities remain independent current directions.
velocity = np.array([1., 1., -2., -2.])
h = 2*np.pi*np.diag(np.abs(velocity))
k = np.diag(np.sign(velocity))/(2*np.pi)
repeated = np.linalg.eigvalsh(np.diag(np.sqrt(np.diag(h)))@k@np.diag(np.sqrt(np.diag(h))))
assert np.allclose(repeated, [-2, -2, 1, 1])
report = dict(records=records, repeated_free_speeds=repeated.tolist(), passed=True)
Path(__file__).with_name('multivelocity_results.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
