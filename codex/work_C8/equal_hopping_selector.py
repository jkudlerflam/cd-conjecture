"""Exact equal-hopping selector. No claims about a large-c payload.

Small complete spectra check the domain-wall mapping, all constraint sectors
check the penalty bound, and independent entropy/Casimir estimates check c=1.
numpy and scipy only. Set OPENBLAS_NUM_THREADS=1 before running.
"""
import itertools
import json
import numpy as np
from scipy.sparse import coo_matrix


def spin_hamiltonian(L, mu=0., a=1., violations=0):
    """Arrow representation b_i with a_i=b_(i-1) xor epsilon_(i-1)."""
    rows, cols, data = [], [], []
    for s in range(1 << L):
        bits = [(s >> j) & 1 for j in range(L)]
        left = [bits[(j-1) % L] ^ ((violations >> ((j-1) % L)) & 1)
                for j in range(L)]
        diag = mu * sum(left[j] == 1 and bits[j] == 0 for j in range(L))
        rows.append(s); cols.append(s); data.append(diag)
        for j in range(L):
            if not ((violations >> j) & 1) and left[j] != bits[(j+1) % L]:
                rows.append(s ^ (1 << j)); cols.append(s); data.append(-a)
    return coo_matrix((data, (rows, cols)), shape=(1 << L,)*2).tocsr()


def wall_spectrum(L, mu=0., a=1.):
    """Even W and both boson twists. JW twist is -g for even W."""
    out = []
    for g in [1, -1]:
        phase = .5 if g == 1 else 0.
        eps = -2*a*np.cos(2*np.pi*(np.arange(L)+phase)/L) + mu/2
        for W in range(0, L+1, 2):
            out.extend(sum(eps[j] for j in c) for c in itertools.combinations(range(L), W))
    return np.sort(out)


def spectrum_checks():
    checks = []
    for L in [5, 6, 7, 8, 9]:
        mu = .37
        H = spin_hamiltonian(L, mu).toarray()
        err = float(np.max(np.abs(np.linalg.eigvalsh(H)-wall_spectrum(L, mu))))
        assert err < 1e-10
        checks.append(dict(L=L, dimension=len(H), mu=mu, spectrum_error=err))
    return checks


def penalty_checks():
    # Test the sharper norm inequality before adding K times the defect count.
    out = []
    for L in [4, 5, 6]:
        mu = .37; bound = 2+abs(mu)
        H0 = spin_hamiltonian(L, mu).toarray()
        e0 = np.linalg.eigvalsh(H0)[0]
        max_ratio = 0.; min_margin = float('inf')
        for defects in range(1, 1 << L):
            k = defects.bit_count()
            Hd = spin_hamiltonian(L, mu, violations=defects).toarray()
            norm = np.max(np.abs(np.linalg.eigvalsh(Hd-H0)))
            max_ratio = max(max_ratio, norm/k)
            min_margin = min(min_margin, np.linalg.eigvalsh(Hd)[0]-e0+bound*k)
        assert max_ratio <= bound+1e-10
        assert min_margin >= -1e-10
        out.append(dict(L=L, checked_patterns=(1 << L)-1,
                        analytic_norm_bound_per_defect=bound,
                        max_norm_per_defect=float(max_ratio),
                        min_energy_lower_bound_margin=float(min_margin)))
    return out


def ground_states(L):
    """Positive XX wall ground state and its legal-arrow isometric image.

    Occupied momenta are consecutive. The absolute Slater determinant is
    the sine Vandermonde product, independent of a common momentum shift.
    L is a multiple of four, W=L/2 is even.
    """
    assert L % 4 == 0
    W = L//2
    xx = np.zeros(1 << L)
    arrow = np.zeros(1 << L)
    mask = (1 << L)-1
    pairs = list(itertools.combinations(range(W), 2))
    combos = itertools.combinations(range(L), W)
    while True:
        batch = list(itertools.islice(combos, 4096))
        if not batch:
            break
        pos = np.asarray(batch, dtype=np.int64)
        amp = np.ones(len(pos))/L**(W/2)
        for i, j in pairs:
            amp *= 2*np.sin(np.pi*(pos[:, j]-pos[:, i])/L)
        wall_ids = np.sum(1 << pos, axis=1)
        xx[wall_ids] = amp
        # Fix b_0=0; w_j=b_(j-1) xor b_j for j=1,...,L-1.
        bit = np.zeros(len(pos), dtype=np.int64)
        spin_ids = np.zeros(len(pos), dtype=np.int64)
        for j in range(1, L):
            bit ^= (wall_ids >> j) & 1
            spin_ids |= bit << j
        arrow[spin_ids] = amp/np.sqrt(2)
        arrow[mask ^ spin_ids] = amp/np.sqrt(2)
    assert abs(xx@xx-1) < 1e-10
    assert abs(arrow@arrow-1) < 1e-10
    return xx, arrow


def entropy_values(singular_values):
    p = singular_values**2
    p = p[p > 1e-15]
    return float(-np.sum(p*np.log(p)))


def xx_entropy(psi, L, ell):
    return entropy_values(np.linalg.svd(psi.reshape(1 << (L-ell), 1 << ell),
                                       compute_uv=False))


def physical_entropy(psi, L, ell):
    """Physical sites 1,...,ell duplicate boundary bits b_0 and b_ell.

    Trace gives four blocks with those endpoint bits fixed. The internal
    indices are b_1,...,b_(ell-1). This avoids a 4^L state vector.
    """
    total = 0.
    left = np.arange(1 << (ell-1), dtype=np.int64) << 1
    right = np.arange(1 << (L-ell-1), dtype=np.int64) << (ell+1)
    for b0, be in itertools.product([0, 1], repeat=2):
        ids = left[None, :] | right[:, None] | b0 | (be << ell)
        total += entropy_values(np.linalg.svd(psi[ids], compute_uv=False))
    return total


def entropy_checks():
    out = []
    for L in [8, 12, 16, 20]:
        xx, psi = ground_states(L)
        W = L//2
        expected = -2*np.sin(np.pi*W/L)/np.sin(np.pi/L)
        if L <= 16:
            residual = float(np.linalg.norm(spin_hamiltonian(L)@psi-expected*psi))
            assert residual < 1e-9
        else:
            residual = None
        ell_values = list(range(max(2, L//4), L//2+1))
        ent = [physical_entropy(psi, L, ell) for ell in ell_values]
        ent_xx = [xx_entropy(xx, L, ell) for ell in ell_values]
        x = np.log((L/np.pi)*np.sin(np.pi*np.array(ell_values)/L))
        fit = np.polyfit(x, ent, 1)
        fit_xx = np.polyfit(x, ent_xx, 1)
        row = dict(L=L, W=W, physical_dimension_per_site=4,
                   ground_vector_residual=residual, ell=ell_values,
                   entropy=ent, xx_control_entropy=ent_xx,
                   c_interval=float(3*fit[0]), c_xx_interval=float(3*fit_xx[0]),
                   half_entropy=ent[-1], xx_half_entropy=ent_xx[-1])
        if out:
            old = out[-1]
            row['c_half_size'] = float(3*(ent[-1]-old['half_entropy'])/np.log(L/old['L']))
            row['c_xx_half_size'] = float(3*(ent_xx[-1]-old['xx_half_entropy'])/np.log(L/old['L']))
        out.append(row)
        print(json.dumps(row), flush=True)
    return out


def casimir_checks():
    sizes = np.array([40, 48, 64, 80, 96, 128, 160], dtype=float)
    e = -2/np.sin(np.pi/sizes)
    # Include the leading analytic correction, no fitted velocity.
    fit = np.polynomial.polynomial.polyfit(1/sizes**2, e/sizes, 2)
    velocity = 2.
    c = -6*fit[1]/(np.pi*velocity)
    first_twist_gap = 2*np.tan(np.pi/(2*sizes))
    return dict(sizes=sizes.tolist(), energies=e.tolist(), velocity=velocity,
                c_casimir=float(c), bulk_energy_fit=float(fit[0]),
                bulk_energy_exact=-2/np.pi,
                twist_gap_times_L=(first_twist_gap*sizes).tolist(),
                twist_scaling_dimension_limit=.25)


def main():
    out = dict(spectra=spectrum_checks(), penalties=penalty_checks(),
               entropy=entropy_checks(), casimir=casimir_checks())
    with open('equal_hopping_selector_checks.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out['casimir'], indent=2))


if __name__ == '__main__':
    main()
