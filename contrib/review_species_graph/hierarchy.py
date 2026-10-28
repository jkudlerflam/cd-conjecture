"""Scale hierarchy and local-dimension accounting for the species-graph candidate.

(1) Explicit p-dependence of every first-order quantity that is fixed by the
    graph and the port encoding.
(2) A conservative LOWER bound on the core alphabet from components that
    are already explicit in contrib/controller_instantiation (circuit record
    D_exact, preloaded dispatcher factor 88*301 = 26,488, fifteen-letter clock,
    three static marks, at least two neutral flag qubits FAIL and READY).
    The program generator, the length counter, the format verifier and the
    readout shuttle are NOT included, so the true C is larger.
(3) The smallest prime p at which (N+1) > d-1 under that lower bound, and
    under the lower bound times extra factors 10^2, 10^4, 10^6.
(4) Asymptotic orders of the gaps that the coupling s must undercut.

Run: python3 hierarchy.py   (instant)
"""
import json, math
from pathlib import Path


def primes(upto):
    s = bytearray([1]) * (upto + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(upto ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(upto + 1) if s[i]]


def graph(p):
    u = p * p + p + 1
    return u, 2 * u, (p + 1) * u


def d_exact(p):
    u, q, N = graph(p)
    ell = math.ceil(math.log2(max(u, N)))
    m = max(1, math.ceil(math.log2(q)), 2 * ell)
    return q + 3 + 2 * (q + 1) * m


COMPONENT_FACTOR = 26488 * 15 * 3 * 4   # dispatcher x clock x marks x (FAIL,READY)


def d_core_lower(p, extra=1.0):
    return COMPONENT_FACTOR * d_exact(p) * extra


def first_violation(extra):
    """Smallest p (searched over a geometric grid of primes and then refined) with
    N+1 > 96 d_core + 48 - 1.  Works with floats for very large p."""
    def ok(p):
        u, q, N = graph(p)
        return N + 1 > 96 * d_core_lower(p, extra) + 47
    lo, hi = 2, 2
    while not ok(hi):
        lo, hi = hi, hi * 2
    while hi - lo > max(1, hi // 10 ** 6):
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    p = hi
    u, q, N = graph(p)
    d = 96 * d_core_lower(p, extra) + 48
    return dict(extra_factor=extra, p_first_approx=p, N_plus_1=float(N + 1), d=float(d),
                C_implied=d_core_lower(p, extra) / (q * (1 + math.log2(q))))


def table():
    rows = []
    for p in [2, 3, 5, 7, 11, 101, 1009, 10007]:
        u, q, N = graph(p)
        eps = 1 / (2 * math.sqrt(u * N))
        eta = eps ** 2 * (1 - 1 / (2 * N))
        dc = d_core_lower(p)
        rows.append(dict(p=p, q=q, N=N, eps=eps, eta=eta,
                         payload_mismatch=math.pi / (6 * (N + 1)),
                         mixed_mismatch=math.pi / 6,
                         D_exact=d_exact(p), d_core_lower=float(dc),
                         C_lower=dc / (q * (1 + math.log2(q))),
                         d_lower=float(96 * dc + 48), c=N + 1,
                         ratio_c_over_d_lower=(N + 1) / (96 * dc + 48)))
    return rows


def main():
    out = dict(component_factor=COMPONENT_FACTOR, table=table(),
               first_violation=[first_violation(x) for x in (1.0, 1e2, 1e4, 1e6)])
    # orders of the gaps (n = core segment length, at least the circuit workspace)
    out['orders'] = {
        'circuit_sites_p2_actual': 81342,
        'cell_length_lower_order': 'n >= 32 A with A = L^2 w, L = N+4q ~ p^3: n >~ 32 p^6 log p; '
                                   'dispatcher layer fits in 2,310,035 A^3 ~ p^18 log^3 p sites',
        'history_length': 'C_n = (2n+4) 2^n + n - 4',
        'clock_gap': '>= J/(4 C_n^2) ~ 4^-n / n^2',
        'fail_acceptance': 'a_n >= eta_F / C_n ; B_n = eta_F n',
        'core_gap': 'Delta_n >= a_n gamma_n/(a_n+B_n+gamma_n) ~ 8^-n / n^4',
        'wrapper_cell_gap': 'tagged-vacancy bound ~ Delta_n / poly(n)',
        'compensator_norm': '||E_L|| <= (1/r1+1/r2)/sin^2(1/(2b)) ~ 8 b^2',
        'normalised_link_norm_B_p': '~ ||E_L|| ||E_R|| 4uN/(p_L p_R) ~ 256 b^4 p^5/(p_L p_R)',
        'coupling_condition': 's << Delta_base * eps_p /(C_p B_p^2): doubly exponentially small in p',
        'payload_RG_condition': 'N g log(1/delta_pay) << 1 with g ~ s B_p^2/Delta_base, i.e. g << 1/(p^3 log p)',
    }
    print(json.dumps(out, indent=1))
    Path(__file__).with_name('hierarchy_results.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
