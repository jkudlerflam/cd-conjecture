"""An initialized qubit helix selects the occupied length of two cells."""
import itertools
import json
import numpy as np
from vacancy_parent import charge_hamiltonian


def plus(theta):
    return np.array([1., np.exp(1j*theta)])/np.sqrt(2)


def tensors(xs):
    out = np.array([1.])
    for x in xs:
        out = np.kron(out, x)
    return out


def filter_data(k):
    phase = 2*np.pi/(k+2)
    singlet = np.array([0, 1, -np.exp(-1j*phase), 0])/np.sqrt(2)
    h = np.outer(singlet, singlet.conj())
    a = np.eye(2)-np.outer(plus(0), plus(0).conj())
    b = np.eye(2)-np.outer(plus((k-1)*phase), plus((k-1)*phase).conj())
    kap = 1/((k+2)*(2*k+3))
    beta = np.sin(np.pi/(k+2))**2
    epsilon = kap*beta/((1+np.sqrt(kap))**2+beta)
    return phase, h, a, b, epsilon


def filter_core(r, h, a, b):
    out = np.kron(a, np.eye(2**(r-1)))+np.kron(np.eye(2**(r-1)), b)
    for j in range(r-1):
        out += np.kron(np.kron(np.eye(2**j), h), np.eye(2**(r-j-2)))
    return out


def occupied_sector(n, r, U, w, scale, h, a, b):
    p = n-r
    if p == 0:
        return scale*filter_core(r, h, a, b), [()], np.array([1.]), 0.
    holes, hc, chi, ec = charge_hamiltonian(r, p, U, w)
    d = 2**r
    out = np.kron(hc, np.eye(d))
    initialized = np.kron(a, np.eye(d//2))+np.kron(np.eye(d//2), b)
    out = out.astype(complex)+scale*np.kron(np.eye(len(holes)), initialized)
    for j in range(r-1):
        indicator = []
        for hs in holes:
            occ = [l for l in range(n) if l not in hs]
            indicator.append(int(occ[j+1] == occ[j]+1))
        payload = np.kron(np.kron(np.eye(2**j), h), np.eye(2**(r-j-2)))
        out += scale*np.kron(np.diag(indicator), payload)
    return out, holes, chi, ec[0]


def one_case(k):
    U, w = 1., .1
    phase, h, a, b, epsilon = filter_data(k)
    scale = 2*U/epsilon
    core_tests = []
    for r in range(2, k+3):
        hr = filter_core(r, h, a, b)
        e = np.linalg.eigvalsh(hr)
        expected = tensors([plus(j*phase) for j in range(r)])
        left_only = hr-np.kron(np.eye(2**(r-1)), b)
        le = np.linalg.eigvalsh(left_only)
        kap_r = 1/(r*(2*r-1))
        beta_r = np.sin((r-k)*phase/2)**2
        individual = kap_r*beta_r/((1+np.sqrt(kap_r))**2+beta_r)
        assert np.linalg.norm(left_only@expected) < 1e-12
        assert le[1] >= kap_r-1e-12
        if r != k:
            assert e[0] >= individual-1e-12 and e[0] >= epsilon-1e-12
        else:
            assert abs(e[0]) < 1e-12 and e[1] > 0
        core_tests.append(dict(r=r, ground_energy=float(e[0]),
                               individual_bound=float(individual), left_gap=float(le[1]),
                               left_gap_bound=kap_r))
    number_tests = []
    for n in (k+1, k+2):
        energies = []
        target_error = 0.
        for r in range(2, n+1):
            hr, holes, chi, ec = occupied_sector(n, r, U, w, scale, h, a, b)
            er = np.linalg.eigvalsh(hr)
            if r == k:
                frame = np.kron(chi, tensors([plus(j*phase) for j in range(r)]))
                target_error = np.linalg.norm(hr@frame-ec*frame)
                assert target_error < 1e-10 and abs(er[0]-ec) < 1e-10
                target_energy = er[0]
            else:
                assert er[0] >= 2*U-w*(n-1)-1e-10
            energies.append(dict(r=r, dimension=hr.shape[0], ground_energy=float(er[0])))
        margin = min(row['ground_energy']-target_energy for row in energies if row['r'] != k)
        assert margin >= U-w*(n-1)-1e-10
        number_tests.append(dict(length=n, sectors=energies, target_code_error=float(target_error),
                                 measured_selection_margin=float(margin),
                                 guaranteed_selection_margin=U-w*(n-1)))
    return dict(k=k, phase=float(phase), U=U, w=w, uniform_filter_bound=float(epsilon),
                filter_scale=float(scale), cores=core_tests, number_selection=number_tests)


if __name__ == '__main__':
    results = [one_case(k) for k in (3, 4)]
    with open('length_filter_checks.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results, indent=2))
