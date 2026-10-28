"""Checks for contrib/minimal_d0 (numpy only, about a minute).
(a) Span of conjugates of the Bell check Pi_Omega (x) 1_c by unitaries on (s,c), at q = 2 and 4.
(b) Exact q = 2 synthesis of the route-B dressed-Bell cell from one static Bell term on (p,s),
    positive time weights on frames SWAP(s,c).Pauli_c, and one non-port term on (s,c).
(c) q = 4: Bell-diagonal (p,s) junk keeps the dressed-Bell code exactly degenerate, with
    compressions proportional to single Majoranas (flavour-dependent alpha)."""
import numpy as np
rng = np.random.default_rng(7)
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex)
P1 = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
def kron(*a):
    o = np.eye(1)
    for x in a: o = np.kron(o, x)
    return o
def haar(n):
    A = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)); Q, R = np.linalg.qr(A)
    return Q*(np.diag(R)/abs(np.diag(R)))
def hvec(H): return np.concatenate([H.real.ravel(), H.imag.ravel()])
def bell(q): v = np.eye(q).reshape(-1)/np.sqrt(q); return np.outer(v, v)

# (a) span ranks
for q, nsamp in [(2, 300), (4, 1400)]:
    D = 2*q*q; menu = np.kron(bell(q), I2)
    vecs = [hvec(np.eye(D))] + [hvec((W := np.kron(np.eye(q), haar(2*q))).conj().T @ menu @ W) for _ in range(nsamp)]
    S = np.linalg.svd(np.array(vecs).T, compute_uv=False)
    print(f"(a) q={q}: rank {int((S > 1e-8*S[0]).sum())}, predicted 4q^4-(4q^2-1)-(q^2-1) = {4*q**4-(4*q*q-1)-(q*q-1)}")

# (b) q = 2 exact synthesis. Order of registers: p, s, c.
SWsc = kron(I2, np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]))
Pi_ps = np.kron(bell(2), I2)
Delta, g = 1.0, 0.37
f1 = 0.05; f2 = f1 + 2*g/Delta; f = (f1 + f2)/2; w0 = 1.0
# Bell frames of (p,c): W = SWAP(s,c) then Pauli on c; Pauli P maps Phi+ to the four Bell states
frames = {'Phi+': ('I', f), 'Phi-': ('Z', f), 'Psi+': ('X', f1), 'Psi-': ('Y', f2)}
H = Delta*w0*(np.eye(8) - Pi_ps).astype(complex)
for name, (P, w) in frames.items():
    W = SWsc @ kron(I2, I2, P1[P])
    H += Delta*w*(np.eye(8) - W.conj().T @ Pi_ps @ W)
H += g*(kron(I2, X, X) - kron(I2, Y, Y))/2          # the non-port (s,c) event
# target: dressed Bell cell H_L = Delta*w0*(1-Pi) + g(B+B^dag) + const
Om = np.eye(2).reshape(-1)/np.sqrt(2)
B = sum(np.kron(np.outer(np.kron(F, I2) @ Om, Om), G) for F, G in [(X, X), (Y, Y)])
HL = Delta*w0*(np.eye(8) - Pi_ps) + g*(B + B.conj().T)
diff = H - HL; const = np.trace(diff).real/8
print("(b) q=2 synthesis residual (after removing a constant):", np.linalg.norm(diff - const*np.eye(8)))
ev, V = np.linalg.eigh(H); G0 = V[:, :2]
print("    lowest levels:", np.round(ev[:3] - ev[0], 12))
Vom = np.kron(Om.reshape(-1, 1), I2)                  # |Omega> (x) psi, logical basis = qubit c
U_, s_, Vh_ = np.linalg.svd(G0.conj().T @ Vom); Iso = G0 @ (U_ @ Vh_)
for nm, F, Gm in [('X_p', kron(X, I2, I2), X), ('Y_p', kron(Y, I2, I2), Y)]:
    C = Iso.conj().T @ F @ Iso; a = np.trace(Gm.conj().T @ C).real/2
    print(f"    compression of {nm} = {a:+.6f} x logical Majorana, residual {np.linalg.norm(C - a*Gm):.2e}")
print("    (Z_p is even under the grading; its compression is the logical parity bilinear, not used by the bond)")

# (c) q = 4: Pauli-string channels odd under Z on the first port qubit, 4 logical modes (JW), Bell-diagonal junk
q = 4; n = 4
def jw(k, which):   # Majorana on n qubits
    ops = [Z]*k + [X if which == 0 else Y] + [I2]*(n-k-1); return kron(*ops)
gam = [jw(k, w) for k in range(n) for w in (0, 1)]
chan = [kron(P1[A], P1[Q]) for A in 'XY' for Q in 'IXYZ']           # r = 8 odd channels
Om4 = np.eye(q).reshape(-1)/np.sqrt(q)
Dl = 2**n
def bellP(Pm):   # projector onto (1 (x) P)|Omega>
    v = np.kron(np.eye(q), Pm) @ Om4; return np.outer(v, v.conj())
junk = sum(rng.uniform(0.0, 0.6)*bellP(kron(P1[a], P1[b])) for a in 'IXYZ' for b in 'IXYZ' if a+b != 'II')
Hps = 1.0*(np.eye(q*q) - bellP(np.eye(q))) + junk
Bm = sum(np.kron(np.outer(np.kron(F, np.eye(q)) @ Om4, Om4), G) for F, G in zip(chan, gam))
H4 = np.kron(Hps, np.eye(Dl)) + 0.6*(Bm + Bm.conj().T)
ev, V = np.linalg.eigh(H4)
print("(c) q=4: lowest 18 levels minus ground:", np.round(ev[:18] - ev[0], 10))
G0 = V[:, :Dl]
# code isometry: map logical psi to ground vector via overlap with |Omega> psi
Vom = np.kron(Om4.reshape(-1, 1), np.eye(Dl))           # |Omega> (x) psi
M = G0.conj().T @ Vom; U_, s_, Vh_ = np.linalg.svd(M); Iso = G0 @ (U_ @ Vh_)   # isometry aligned with logical basis
worst = 0
for F, G in zip(chan, gam):
    C = Iso.conj().T @ np.kron(np.kron(F, np.eye(q)), np.eye(Dl)) @ Iso
    a = np.trace(G.conj().T @ C).real/Dl
    worst = max(worst, np.linalg.norm(C - a*G))
    print(f"    alpha = {a:+.6f}", end="")
print(f"\n    max residual of compression minus alpha*Majorana: {worst:.2e}")
