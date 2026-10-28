"""Re-implementation of the 'exactly 2^N initialized inputs' mechanism.

Part A.  Query signatures.  The verifier of EXISTENTIAL_CORE reads records
only through a TYPE query and indexed species-bit queries.  We build the
complete record alphabet D = q+3+2(q+1)m, compute every record's full query
signature, and count the accepted words of a small register under (i) the
original format (a neutral slot required to hold an ordinary bit 0) and
(ii) the repaired format (neutral slots must hold the unique blank).

Part B.  A small Feynman-Kitaev history Hamiltonian (abstract clock, not
translation invariant) for: format check -> copy violation to FAIL ->
encode -> readout -> decode -> uncheck.  It compares an all-time FAIL
penalty with a terminal-only one, checks the kernel dimension, the
acceptance-gap lemma gap >= a*gamma/(a+B+gamma), and the readout
compression p_read * V^dag F V on the kernel.

numpy only.  Run: python3 input_selection.py  (seconds)
"""
import itertools, json
from pathlib import Path
import numpy as np

RES = {}


# ------------------------------------------------------------------ Part A
def alphabet(q, m):
    return ([('D', s) for s in range(q)] + [('blank',), ('bit', 0), ('bit', 1)]
            + [('H', j, c) for j in range(m) for c in (0, 1)]
            + [('T', j, c, s) for j in range(m) for c in (0, 1) for s in range(q)])


def signature(rec, nbits):
    typ = {'D': 0, 'blank': 1, 'bit': 2, 'H': 3, 'T': 4}[rec[0]]
    bits = tuple(((rec[1] >> j) & 1) if rec[0] == 'D' else 0 for j in range(nbits))
    return (typ,) + bits


def part_a():
    q, m = 3, 2
    nb = int(np.ceil(np.log2(q)))
    A = alphabet(q, m)
    sig = {r: signature(r, nb) for r in A}
    classes = {}
    for r, s in sig.items():
        classes.setdefault(s, []).append(r)
    aliased = [c for c in classes.values() if len(c) > 1]
    # register: two occupation slots (allowed {D_0,D_1} and {D_1,D_2}), one charged
    # reference slot (D_2), one neutral slot.
    occ_allowed = [{('D', 0), ('D', 1)}, {('D', 1), ('D', 2)}]
    ref = ('D', 2)

    def accept(word, neutral_target):
        ok = all(sig[word[i]] in {sig[x] for x in occ_allowed[i]} for i in range(2))
        ok &= sig[word[2]] == sig[ref]
        ok &= sig[word[3]] == sig[neutral_target]
        return ok

    counts = {}
    for label, tgt in [('original_neutral_bit0', ('bit', 0)), ('repaired_blank', ('blank',))]:
        counts[label] = sum(accept(w, tgt) for w in itertools.product(A, repeat=4))
    RES['query_signatures'] = dict(q=q, m=m, alphabet=len(A), signature_classes=len(classes),
                                   aliased_classes=[[str(x) for x in c] for c in aliased],
                                   accepted_words=counts, intended=4)
    print(json.dumps(RES['query_signatures'], indent=1))
    assert counts['repaired_blank'] == 4 and counts['original_neutral_bit0'] == 8


# ------------------------------------------------------------------ Part B
def part_b(eta=1.0, J=1.0):
    # registers: two data qutrits (valid values 0,1; value 2 malformed), FAIL bit,
    # two violation scratch bits v1,v2, READY bit.  Basis index (d1,d2,F,v1,v2,R).
    dims = [3, 3, 2, 2, 2, 2]
    D = int(np.prod(dims))
    idx = lambda t: int(np.ravel_multi_index(t, dims))
    states = list(itertools.product(*[range(k) for k in dims]))

    def perm(f):
        U = np.zeros((D, D))
        for s in states:
            U[idx(f(s)), idx(s)] = 1
        return U

    def X_v(i):   # v_i ^= [d_i == 2]
        return perm(lambda s: tuple(s[:3]) + ((s[3] ^ (s[0] == 2)) if i == 0 else s[3],
                                              (s[4] ^ (s[1] == 2)) if i == 1 else s[4], s[5]))
    OR_to_F = perm(lambda s: (s[0], s[1], s[2] ^ int(s[3] or s[4]), s[3], s[4], s[5]))
    READY = perm(lambda s: s[:5] + (s[5] ^ 1,))
    # encoding: a fixed random unitary on the valid 4-dim data space, identity elsewhere
    rng = np.random.default_rng(7)
    Z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    Uv, _ = np.linalg.qr(Z)
    Uenc = np.eye(D, dtype=complex)
    valid = [(a, b) for a in (0, 1) for b in (0, 1)]
    for rest in itertools.product(range(2), range(2), range(2), range(2)):
        ids = [idx((a, b) + rest) for a, b in valid]
        Uenc[np.ix_(ids, ids)] = Uv
    I = np.eye(D)
    gates = [X_v(0), X_v(1), OR_to_F, X_v(0), X_v(1),      # check, keep FAIL
             Uenc, READY, I, I, READY,                        # encode, readout window
             Uenc.conj().T, X_v(0), X_v(1), OR_to_F, X_v(0), X_v(1)]  # decode, uncheck
    T = len(gates)
    C = T + 1
    W = [I.astype(complex)]
    for g in gates:
        W.append(g @ W[-1])
    proj = lambda reg, val: np.diag([1.0 if s[reg] == val else 0.0 for s in states])
    clock = lambda t1, t2: np.outer(np.eye(C)[t1], np.eye(C)[t2])
    Hprop = np.zeros((C * D, C * D), complex)
    for t, g in enumerate(gates):
        Hprop += 0.5 * J * (np.kron(clock(t, t), I) + np.kron(clock(t + 1, t + 1), I)
                            - np.kron(clock(t + 1, t), g) - np.kron(clock(t, t + 1), g.conj().T))
    Hinit = np.kron(clock(0, 0), proj(2, 1) + proj(3, 1) + proj(4, 1) + proj(5, 1))
    Hbase = Hprop + Hinit
    A_all = eta * sum(np.kron(clock(t, t), proj(2, 1)) for t in range(C))
    A_term = eta * np.kron(clock(T, T), proj(2, 1))
    out = {}
    ev0 = np.linalg.eigvalsh(Hbase)
    gamma = float(ev0[ev0 > 1e-9][0])
    out['initialized_kernel_dim'] = int((ev0 < 1e-9).sum())   # = 9 data inputs
    for name, A in [('all_time', A_all), ('terminal_only', A_term)]:
        ev, vec = np.linalg.eigh(Hbase + A)
        k = int((ev < 1e-9).sum())
        out[name] = dict(kernel_dim=k, gap=float(ev[k]) if k < len(ev) else None)
    # acceptance lemma: a >= eta * (#FAIL-display times)/C ; B = ||A||
    invalid = [idx((d1, d2, 0, 0, 0, 0)) for d1 in range(3) for d2 in range(3) if 2 in (d1, d2)]
    PF = proj(2, 1)
    nfail = min(sum(np.linalg.norm(PF @ W[t][:, x]) ** 2 for t in range(C)) for x in invalid)
    a = eta * nfail / C; B = eta; bound = a * gamma / (a + B + gamma)
    out['acceptance_bound'] = dict(gamma=gamma, a=a, B=B, bound=bound,
                                   holds=bool(out['all_time']['gap'] >= bound - 1e-12))
    # readout compression on the accepted kernel: Pi_READY (x) F, F a matrix unit on d1
    ev, vec = np.linalg.eigh(Hbase + A_all)
    K = vec[:, ev < 1e-9]
    F = perm(lambda s: ((s[0] + 1) % 2 if s[0] < 2 else s[0],) + s[1:])   # flip d1 within {0,1}
    F_read = np.kron(np.eye(C), proj(5, 1) @ F)
    comp = K.conj().T @ F_read @ K
    # predicted: (|S|/C) * V^dag F V on the valid inputs, with V = Uenc
    hist = lambda x: sum(np.kron(np.eye(C)[t], W[t] @ np.eye(D)[x]) for t in range(C)) / np.sqrt(C)
    Hs = np.array([hist(idx((a_, b_, 0, 0, 0, 0))) for a_, b_ in valid]).T
    S = [t for t in range(C) if abs((proj(5, 1) @ W[t])[idx((0, 0, 0, 0, 0, 1)), idx((0, 0, 0, 0, 0, 0))]) > 0.5]
    pred = (len(S) / C) * (Uv.conj().T @ np.kron(np.array([[0, 1], [1, 0]]), np.eye(2)) @ Uv)
    got = Hs.conj().T @ F_read @ Hs
    out['readout'] = dict(active_times=len(S), C=C, err=float(np.abs(got - pred).max()),
                          kernel_spans_valid_histories=float(np.linalg.norm(K @ (K.conj().T @ Hs) - Hs)))
    RES['tiny_history'] = out
    print(json.dumps(out, indent=1))
    assert out['initialized_kernel_dim'] == 9
    assert out['all_time']['kernel_dim'] == 4 and out['terminal_only']['kernel_dim'] == 9
    assert out['acceptance_bound']['holds'] and out['readout']['err'] < 1e-12


if __name__ == '__main__':
    part_a(); part_b()
    Path(__file__).with_name('input_selection_results.json').write_text(json.dumps(RES, indent=1) + '\n')
