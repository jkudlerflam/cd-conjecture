"""Closed-book C10 verification. Only reads/writes files created for this task."""
import json
import itertools
import numpy as np
from scipy.linalg import block_diag, eigh

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1]).astype(complex)

def kron(*xs):
    out = np.ones((1, 1), complex)
    for x in xs:
        out = np.kron(out, x)
    return out

def norm(x):
    return float(np.linalg.norm(x, 2))

def small(x):
    return float(np.real_if_close(x))

m = 4
positions = [2, 3, 2]
n = len(positions)
letters = ['M'] + [f'A{j}' for j in range(1, m+1)] + [f'R{j}' for j in range(1, m+1)] + [f'H{t}' for t in range(n)]
index = {x: j for j, x in enumerate(letters)}
allowed = set()
for j in range(1, m):
    allowed.add((f'A{j}', f'A{j+1}'))
    allowed.add((f'R{j}', f'R{j+1}'))
allowed.update([('M', 'A1'), (f'R{m}', 'M')])
for t, j in enumerate(positions):
    allowed.add(('M' if j == 1 else f'A{j-1}', f'H{t}'))
    allowed.add((f'H{t}', 'M' if j == m else f'R{j+1}'))

rules = []
for t, (j, jp) in enumerate(zip(positions[:-1], positions[1:])):
    if jp == j+1:
        rules.append(((f'H{t}', f'R{j+1}'), (f'A{j}', f'H{t+1}')))
    elif jp == j-1:
        rules.append(((f'A{j-1}', f'H{t}'), (f'H{t+1}', f'R{j}')))
    else:
        raise ValueError('Head must move one site per transition.')

legal_words = []
violation_histogram = {}
legal_edges = set()
# Include interior M labels, which violate selector-track synchronization.
for word in itertools.product(letters, repeat=m):
    seq = ('M',) + word + ('M',)
    violations = sum((seq[j], seq[j+1]) not in allowed for j in range(m+1)) + word.count('M')
    violation_histogram[violations] = violation_histogram.get(violations, 0) + 1
    if violations == 0:
        legal_words.append(word)

word_times = {}
for t, j in enumerate(positions):
    word = tuple([f'A{k}' for k in range(1, j)] + [f'H{t}'] + [f'R{k}' for k in range(j+1, m+1)])
    word_times[word] = t
assert set(legal_words) == set(word_times)
for word in legal_words:
    for site in range(m-1):
        pair = word[site:site+2]
        for a, b in rules:
            for old, new in ((a, b), (b, a)):
                if pair == old:
                    target = word[:site] + new + word[site+2:]
                    assert target in word_times
                    legal_edges.add(tuple(sorted((word_times[word], word_times[target]))))
assert legal_edges == {(0, 1), (1, 2)}

# Transfer-matrix count of legal segments of different lengths, without bulk enumeration.
nonmarker = [x for x in letters if x != 'M']
legal_length_counts = {}
for length in range(1, 9):
    counts = {x: int(('M', x) in allowed) for x in nonmarker}
    for _ in range(1, length):
        counts = {y: sum(counts[x] for x in nonmarker if (x, y) in allowed) for y in nonmarker}
    legal_length_counts[length] = sum(v for x, v in counts.items() if (x, 'M') in allowed)
assert legal_length_counts[4] == n
assert all(count == 0 for length, count in legal_length_counts.items() if length != 4)

# Physical data order: left port p, partner s, logical l, idle right port r.
omega = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
PO = np.outer(omega, omega.conj())
P = kron(PO, I, I)
A = np.kron(omega[:, None], np.eye(4, dtype=complex))
Fp = [kron(X, I, I, I), kron(Y, I, I, I)]
logical = [kron(I, I, X, I), kron(I, I, Y, I)]
B = sum(F @ P @ gamma for F, gamma in zip(Fp, logical))
Delta, g = 1.0, 0.2
HL = Delta*(np.eye(16)-P) + g*(B+B.conj().T)
HL_pauli = Delta*(np.eye(16)-P) + (g/2)*(kron(X,I,X,I)+kron(Y,I,Y,I)+kron(I,X,X,I)-kron(I,Y,Y,I))
assert norm(HL-HL_pauli) < 1e-12

SW = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], complex)
W = kron(I, SW, I)
U = [np.eye(16), W, W]
Ud = block_diag(*U)
V = np.vstack(U)/np.sqrt(n)
Lclock = np.array([[1,-1,0],[-1,2,-1],[0,-1,1]], float)/2
D0 = n*Delta*(np.eye(16)-P) + (n*g/2)*(kron(I,X,X,I)-kron(I,Y,Y,I))
D1 = np.zeros((16,16), complex)
D2 = (n*g/2)*(kron(X,X,I,I)+kron(Y,Y,I,I))
D = block_diag(D0, D1, D2)
compressed = V.conj().T @ D @ V
assert norm(compressed-HL) < 1e-12

r = 2
th = .5*np.arctan2(2*g*np.sqrt(r), Delta)
alpha = -np.sin(2*th)/np.sqrt(r)
Eminus = (Delta-np.sqrt(Delta**2+4*g*g*r))/2
G = np.cos(th)*A - (np.sin(th)/np.sqrt(r))*B@A
Ghist = V@G
assert norm(G.conj().T@G-np.eye(4)) < 1e-12
assert norm(HL@G-Eminus*G) < 1e-12
Fclock = [kron(np.eye(n), F) for F in Fp]
ideal_comp_errors = [norm(Ghist.conj().T@Fc@Ghist-alpha*kron(pauli,I)) for Fc, pauli in zip(Fclock,[X,Y])]
history_invariance = [norm(Fc@V-V@F) for Fc,F in zip(Fclock,Fp)]
ideal_residual = norm(D@Ghist-Eminus*Ghist)
# Projection of the residual outside all propagation histories.
history_leakage_of_D = norm((np.eye(48)-V@V.conj().T)@D@Ghist)
Lplus = np.linalg.pinv(Lclock)
second_order = -Ghist.conj().T@D@Ud@kron(Lplus,np.eye(16))@Ud.conj().T@D@Ghist

scans = []
for J in [10., 100., 1000., 10000.]:
    HP = J*Ud@kron(Lclock,np.eye(16))@Ud.conj().T
    vals, vecs = eigh(HP+D)
    Evec = vecs[:,:4]
    overlap = Ghist.conj().T@Evec
    aa, ss, bb = np.linalg.svd(overlap)
    align = bb.conj().T@aa.conj().T
    Ealign = Evec@align
    errors = [norm(Ealign.conj().T@Fc@Ealign-alpha*kron(pauli,I)) for Fc,pauli in zip(Fclock,[X,Y])]
    compressed_ports = [Ealign.conj().T@Fc@Ealign for Fc in Fclock]
    v = norm(D)
    gamma_clock = J*(1-np.cos(np.pi/n))
    scans.append(dict(J=J, lowest_eigenvalues=[float(x) for x in vals[:8]],
        four_state_splitting=float(vals[3]-vals[0]), gap_above_four=float(vals[4]-vals[3]),
        leakage_from_full_history=norm((np.eye(48)-V@V.conj().T)@Evec),
        distance_from_ideal_code=float(np.sqrt(max(0,1-min(ss)**2))),
        port_compression_errors=errors,
        compressed_port_spectra=[[float(x) for x in np.linalg.eigvalsh(F)] for F in compressed_ports],
        port_square_errors=[norm(F@F-alpha**2*np.eye(4)) for F in compressed_ports],
        ground_energy_error=float(vals[0]-Eminus),
        self_energy_bound=(v*v/(gamma_clock-2*v) if gamma_clock>2*v else None)))

# Exact frustration-free Bell preparation, kept as ENERGY, never a legal-pair constraint.
Jff = 7.
HPff = Jff*Ud@kron(Lclock,np.eye(16))@Ud.conj().T
DBell = block_diag(n*Delta*(np.eye(16)-P), np.zeros((16,16)), np.zeros((16,16)))
ffvals, ffvecs = eigh(HPff+DBell)
Gff = V@A
# The Bell projector gives the same rank-one-clock bound as an initial-input check.
init_bound = (n*Delta*Jff)/(2*n*(Jff+(n-1)*n*Delta))
clock_gap = Jff*(1-np.cos(np.pi/n))
ff_report = dict(ground_dimension=int(np.sum(np.abs(ffvals)<1e-9)),
    zero_energy_residual=norm((HPff+DBell)@Gff),
    ground_projector_error=norm(ffvecs[:,:4]@ffvecs[:,:4].conj().T-Gff@Gff.conj().T),
    gap=float(ffvals[4]), lower_bound=float(min(init_bound,clock_gap)))

# Explicit illegal no-head block: R1 R2 R3 R4 has one forbidden left boundary,
# no propagation rules or time energies, hence kappa I on the 16-dimensional data block.
kappa = 20.
nohead = tuple(f'R{j}' for j in range(1,m+1))
seq = ('M',)+nohead+('M',)
assert sum((seq[j],seq[j+1]) not in allowed for j in range(m+1)) == 1
assert not any(nohead[s:s+2] in (a,b) for s in range(m-1) for a,b in rules)
nohead_spectrum = np.linalg.eigvalsh(kappa*np.eye(16))
C = n*abs(g)
shifted_blocks = [D0+C*np.eye(16), D1+C*np.eye(16), D2+C*np.eye(16)]
shifted_minima = [float(np.linalg.eigvalsh(d)[0]) for d in shifted_blocks]
assert min(shifted_minima) > -1e-12

# A pure clock sanity check, including the entire tower.
clockvals = np.linalg.eigvalsh(Jff*Lclock)
clockexpected = Jff*(1-np.cos(np.arange(n)*np.pi/n))
assert np.max(np.abs(clockvals-clockexpected)) < 1e-12

report = dict(
    alphabet=letters, k_c=len(letters), q=2, d0=4*len(letters)*2,
    allowed_pairs=sorted([list(x) for x in allowed]),
    transition_patterns=[[list(a),list(b)] for a,b in rules],
    all_interior_track_strings=len(letters)**m,
    violation_histogram=violation_histogram,
    legal_words=[list(w) for w in legal_words],
    legal_length_counts=legal_length_counts,
    legal_edges=sorted([list(x) for x in legal_edges]),
    route_b=dict(Delta=Delta,g=g,Eminus=Eminus,alpha=alpha,
        pauli_expansion_error=norm(HL-HL_pauli), compressed_energy_error=norm(compressed-HL),
        ideal_port_compression_errors=ideal_comp_errors,
        full_history_port_invariance_errors=history_invariance,
        uniform_history_eigenvector_residual=ideal_residual,
        D_leakage_out_of_history=history_leakage_of_D,
        second_order_energy_coefficients=[float(x) for x in np.linalg.eigvalsh(second_order)],
        clock_energy_norm=norm(D), finite_J=scans),
    frustration_free_Bell_energy=ff_report,
    clock_tower=dict(observed=[float(x) for x in clockvals],expected=[float(x) for x in clockexpected]),
    illegal_sector=dict(kappa=kappa,C=C,shifted_time_block_minima=shifted_minima,
        explicit_nohead_block_min=float(nohead_spectrum[0]),
        all_illegal_lower_bound_correct_marker_labels=kappa,
        all_illegal_lower_bound_shared_boundary_penalties=kappa/2,
        after_marker_counterterm_bound=kappa/2-C),
    scope='Full 48-dimensional legal-sector ED; exhaustive 12^4 classical legality check; explicit 16-dimensional illegal block ED. No full quantum illegal-sector ED. Right port is a free spectator; only the left Route B gadget is implemented.'
)
with open('verification.json','w') as f:
    json.dump(report,f,indent=2)
print(json.dumps(report,indent=2))
