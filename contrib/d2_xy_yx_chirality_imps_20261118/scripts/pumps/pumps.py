#!/usr/bin/env python3
"""Periodic uniform MPS on an oriented ring. Run python3 pumps.py --help.
Written by OpenAI Codex (gpt-6.2-astra, xhigh) to a spec by Claude Opus 5.5, 2026-11-24; self-checks in selfcheck.txt.

The public Ring methods act on complex tensors indexed (physical,left,right).
Tangent columns use C order. Energies are for the whole ring. Translation
moves a spin to its right. Phi_m carries exp(-2 pi i m j/L) on the ket
defect and Ht_n carries exp(+2 pi i n (j+1/2)/L) on the oriented bond.

Ring.tangent(m,E0) returns the dense norm and centered Hamiltonian.
excitations returns (omega, B_columns, diagnostics). ks_data consumes
these columns at m = -n mod L. No spin quantum number is imposed.

A gap of k unchanged sites has transfer E^k. Its one-bond counterpart
F_k is the sum of E^i V E^(k-2-i). Appending a site gives
F_(k+1) = F_k E + E^(k-1) V, with F_0 = F_1 = 0.
Translation fixes the bra defect. Each remaining ket position contributes
two arcs, including their endpoint bonds, and one separate coincident
case closes the ring. Each momentum costs O(L d^2 D^6) for the matrices
and O((d D^2)^3) for diagonalization. Short rings use bounded dyadic
checkpoints. On long rings finite polynomial coefficient extraction
keeps the same asymptotic runtime with O(d^2 D^4) workspace.

The ED oracle is independent of these contractions. It embeds the
ordered two-site matrix by Kronecker products and Fourier transforms
translation orbits. check=1 also compares against literal trace vectors.
"""
import os
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.optimize as opt

if os.environ.get("PUMPS_BACKEND", "numpy").lower() == "cupy":
    import cupy as xp
elif os.environ.get("PUMPS_BACKEND", "numpy").lower() == "numpy":
    xp = np
else:
    raise ValueError("PUMPS_BACKEND must be numpy or cupy")

TAU = 2 * np.pi


def host(x):
    return np.asarray(x) if xp is np else xp.asnumpy(x)


def sync():
    if xp is not np:
        xp.cuda.get_current_stream().synchronize()


def herm(x):
    return (x + x.conj().T) * 0.5


def model(name, delta=0.1875, h=1.783, Jz=0.5, file=None):
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1])
    k = np.kron
    if name == "chir":
        a = (1 + delta) / 2
        op = a*k(X,Y) + (1-a)*k(Y,X) + k(Z,Z) - h/2*(k(Z,I)+k(I,Z))
    elif name == "xxz":
        op = k(X,X) + k(Y,Y) + Jz*k(Z,Z)
    elif name == "tfim":
        op = -k(X,X) - (k(Z,I)+k(I,Z))/2
    elif name == "file":
        op = np.load(file, allow_pickle=False)
    else:
        raise ValueError("unknown model " + name)
    if op.ndim!=2 or len(op)==0 or not np.all(np.isfinite(op)):
        raise ValueError("h2 must be a finite, nonempty matrix")
    d = math.isqrt(op.shape[0])
    if op.shape != (d*d, d*d) or not np.allclose(op, op.conj().T, atol=1e-13, rtol=1e-13):
        raise ValueError("h2 must be a Hermitian d² by d² matrix")
    return xp.asarray(op, dtype=xp.complex128)


def products(a, b):
    return (a[:, None] @ b[None, :]).reshape(-1, a.shape[1], a.shape[2])


def transfer(bra, ket):
    """Closed physical block. Row (bra left, ket left), column (bra right, ket right)."""
    d, D, _ = bra.shape
    out = bra.reshape(d,-1).conj().T @ ket.reshape(d,-1)
    return out.reshape(D,D,D,D).transpose(0,2,1,3).reshape(D*D,D*D)


def apply_op(op, block):
    return (op @ block.reshape(len(block), -1)).reshape(block.shape)


def power_insert(E, V, t, z=1):
    """Return E^t and sum_{j<t} z^j E^j V E^(t-1-j), via segment doubling."""
    if t == 0:
        return xp.eye(len(E), dtype=E.dtype), xp.zeros_like(E)
    P, G = E, V
    k = 1
    for bit in bin(t)[3:]:
        G = G @ P + z**k * (P @ G)
        P = P @ P
        k *= 2
        if bit == "1":
            G = G @ E + z**k * (P @ V)
            P = P @ E
            k += 1
    return P, G


def reverse_powers(E, highest):
    """Stream E^highest,...,I with dyadic checkpoints and O(highest) GEMMs."""
    jumps = [E]
    while 2**len(jumps) <= highest:
        jumps.append(jumps[-1] @ jumps[-1])

    def descend(start, count):
        if count == 1:
            yield start
        else:
            b = (count-1).bit_length()-1
            cut = 1 << b
            yield from descend(start @ jumps[b], count-cut)
            yield from descend(start, cut)
    yield from descend(xp.eye(len(E), dtype=E.dtype), highest+1)


def gap_terms(E,V,L,z,metric_only,polynomial):
    """Yield a bilinear expansion for two complementary gaps.

    Each row gives coefficient, opposite gap, forward gap, its bond sum,
    its one-site-shorter power and whether to include an adjacent contact.
    The direct recursion uses at most seven dyadic levels for L <= 128.
    The long-ring branch has constant workspace even as L tends to infinity.
    """
    K = L-2
    I = xp.eye(len(E),dtype=E.dtype)
    if not polynomial or K<2:
        prev = X = I
        F = xp.zeros_like(E)
        for r,Y in enumerate(reverse_powers(E,K),1):
            yield L*z**r,Y,X,F,prev,r==1,r>1
            if not metric_only and r>1:
                F = F @ E + prev @ V
            prev,X = X,X @ E
        return

    # For any bilinear join J, the coefficient of u^K in
    # J(sum_k (zu)^k X_k, sum_j u^j E^j) pairs k+j=K. Both polynomials
    # are truncated at K. The (2K+1)-point DFT is alias free and exact
    # in exact arithmetic. No transfer eigendecomposition is involved.
    if metric_only:
        Pm = xp.linalg.matrix_power(E,K-1)
    else:
        Pm,FK = power_insert(E,V,K-1)
        tail = E @ FK
        del FK
    P = Pm @ E
    Pplus = P @ E
    scale = float(xp.linalg.norm(P))
    rho = np.exp(-(1+np.log(max(scale,1e-300)))/K)
    count = 2*K+1
    for j in range(count):
        u = rho*np.exp(1j*TAU*j/count)
        a = z*u
        Y = xp.linalg.solve(I-u*E,I-u**(K+1)*Pplus)
        divisor = I-a*E
        X = xp.linalg.solve(divisor,I-a**(K+1)*Pplus)
        F = prev = None
        if not metric_only:
            smaller = X-a**K*P
            # (I-aE) sum a^k F_k = a² V R_{K-2} - a^(K+1) E F_K.
            F = xp.linalg.solve(divisor, a*a*(V @ (smaller-a**(K-1)*Pm))-a**(K+1)*tail)
            prev = a*smaller
        yield L*z/(count*u**K),Y,X,F,prev,True,True


class OpenArc:
    """Virtual-index operations for two open tensors on a cyclic contraction.

    A decorated arc may contain any closed operator insertions. Its free
    endpoint maps and the final joining contraction are independent of the
    number of insertions in either gap.
    """
    def __init__(self,A):
        self.A = A
        self.D = A.shape[1]

    def tensor(self,x):
        return x.reshape((self.D,)*4)

    def left(self,s,x):
        D = self.D
        return (self.A[s] @ x.reshape(D,D,-1)).reshape(x.shape)

    def right(self,x,t):
        v = self.tensor(x).transpose(0,1,3,2) @ self.A[t].conj()
        return v.transpose(0,1,3,2).reshape(x.shape)

    def opposite(self,Y):
        return self.tensor(Y).transpose(3,0,2,1).reshape(Y.shape)

    def join(self,x,y):
        v = self.tensor(x).transpose(0,3,1,2).reshape(x.shape) @ y
        return self.tensor(v).transpose(2,0,1,3).reshape(x.shape)

    def diagonal(self,x):
        return self.tensor(x).transpose(2,0,3,1).reshape(x.shape)


class Ring:
    def __init__(self, A, h2, L):
        self.A = A = xp.asarray(A, dtype=xp.complex128)
        self.h = xp.asarray(h2, dtype=xp.complex128)
        self.L = int(L)
        self.d, self.D, right = A.shape
        if self.L < 2 or self.D != right or self.h.shape != (self.d**2,)*2:
            raise ValueError("incompatible ring or tensor dimensions (L >= 2)")
        self.E = transfer(A,A)
        self.AA = products(A,A)
        self.hAA = apply_op(self.h, self.AA)
        self.V = transfer(self.AA, self.hAA)

    def phase(self, m):
        residue = m % self.L
        if 2*residue == self.L:
            return -1.+0j
        if 2*residue > self.L:
            residue -= self.L
        return np.exp(1j*TAU*residue/self.L)

    def ground_env(self, z=1, centered=0.):
        """One open bra tensor and either the identity or sum z^j h_j.

        Return unnormalized norm, energy and both one-site covectors. The
        derivative of the real total energy uses 2 L times the connected one.
        """
        A, D, d, L = self.A, self.D, self.d, self.L
        C = self.hAA
        V = self.V
        if centered:
            C = apply_op(self.h-centered*xp.eye(d*d),self.AA)
            V = transfer(self.AA,C)
        P, G = power_insert(self.E, V, L-2, z)
        P4 = P.reshape((D,)*4)
        E1 = (P @ self.E).reshape((D,)*4)
        norm_vec = xp.einsum('jbif,sfb->sij', E1, A, optimize=True)
        hv = z*xp.einsum('jbif,sfb->sij', G.reshape((D,)*4), A, optimize=True)
        C = C.reshape(d,d,D,D)
        # Contract the ket block with the empty environment once for both ends.
        edge = xp.einsum('cbif,svfb->svci', P4, C, optimize=True)
        hv += xp.einsum('svci,vjc->sij', edge, A.conj(), optimize=True)
        hv += z**(L-1)*xp.einsum('vsjc,vci->sij', edge, A.conj(), optimize=True)
        norm = float(xp.vdot(A,norm_vec).real)
        energy = L*float(xp.einsum('ij,ji', self.V, P).real)/norm
        return norm, energy, norm_vec, hv

    def norm(self):
        return float(xp.trace(xp.linalg.matrix_power(self.E,self.L)).real)

    def normalized(self):
        # Track the logarithm of the scale for gauges or loads with a large
        # spectral radius. The ordinary path avoids extra reductions on GPU.
        bound = float(xp.linalg.norm(self.E, ord=np.inf))
        if not np.isfinite(bound) or bound == 0:
            raise FloatingPointError("zero or nonfinite transfer matrix")
        if abs(self.L*np.log(bound)) < 500:
            norm = self.norm()
            if np.isfinite(norm) and norm>1e-250:
                return self.A*np.exp(-np.log(norm)/(2*self.L))
        base = self.E/bound
        P = base
        logscale = np.log(bound)
        for bit in bin(self.L)[3:]:
            P = P @ P
            logscale *= 2
            if bit=="1":
                P = P @ base
                logscale += np.log(bound)
            size = float(xp.max(xp.abs(P)))
            if size==0:
                raise FloatingPointError("periodic state has zero norm")
            P /= size
            logscale += np.log(size)
        tr = float(xp.trace(P).real)
        if tr<=0:
            raise FloatingPointError("periodic norm underresolved in this gauge")
        return self.A*np.exp(-(logscale+np.log(tr))/(2*self.L))

    def ks_vector(self, n):
        z = self.phase(n)
        shift = float(xp.trace(self.h).real)/len(self.h) if n%self.L else 0.
        norm, _, _, hv = self.ground_env(z,shift)
        return (self.L**2/TAU * np.exp(1j*np.pi*n/self.L) / np.sqrt(norm) * hv).ravel()

    def weight(self, n, centered=None):
        """Squared norm of the Fourier sum on normalized psi, before L/(2pi).

        centered subtracts a scalar per bond, useful for the energy variance.
        """
        A, d, L = self.A, self.d, self.L
        z = self.phase(n)
        if centered is None:
            centered = float(xp.trace(self.h).real)/len(self.h) if n%L else 0.
        h = self.h - centered*xp.eye(d*d)
        C = apply_op(h, self.AA)
        V = transfer(self.AA,C)
        P = xp.linalg.matrix_power(self.E, L-2)
        corr = xp.einsum('ij,ji', transfer(C,C), P)
        if L == 2:
            rev = apply_op(h.reshape(d,d,d,d).transpose(1,0,3,2).reshape(d*d,d*d), self.AA)
            corr += z*xp.trace(transfer(C,rev))
        else:
            P3, G = power_insert(self.E, V, L-3, z)
            adj = xp.einsum('ij,ji', transfer(products(C,A), products(A,C)), P3)
            corr += 2*(z*adj).real
            corr += z**2*xp.einsum('ij,ji', transfer(C,self.AA), G)
        return L*float(corr.real)/self.norm()

    def tangent(self, m, reference=0., metric_only=False, polynomial=None):
        """Dense norm and centered Hamiltonian. Fix the open bra at site zero.

        For r>0 contract bonds on the clockwise arc [0,r]; exchanging layers
        supplies the complementary arc by Hermitian conjugation. Gap transfers
        have one insertion accumulated by a first-order product recurrence.
        """
        A, E, L, d, D = self.A, self.E, self.L, self.d, self.D
        z = self.phase(-m)
        # Center each bond before contraction to reduce cancellation.
        V = None
        arc = OpenArc(A)
        if polynomial is None:
            polynomial = L>128

        # Block ordering (s,t,ij,kl) makes accumulation contiguous.
        N = xp.zeros((d,d,D*D,D*D), dtype=A.dtype)
        H = None if metric_only else xp.zeros_like(N)
        if not metric_only:
            h = self.h - (reference/L)*xp.eye(d*d)
            C = apply_op(h, self.AA).reshape(d,d,D,D)
            V = self.V if reference==0 else transfer(self.AA,C.reshape(d*d,D,D))
            Pleft = [transfer(A,C[s]) for s in range(d)]
            Qright = [transfer(C[:,t],A) for t in range(d)]
            hh = h.reshape((d,)*4)

            if polynomial:
                contacts = [[transfer(A,apply_op(hh[s,:,:,t],A)) for t in range(d)] for s in range(d)]

        for phase,Y,X,F,prev,contact,separated in gap_terms(E,V,L,z,metric_only,polynomial):
            y = arc.opposite(Y)
            for s in range(d):
                KX = arc.left(s,X)
                if not metric_only and separated:
                    KF = arc.left(s,F) + Pleft[s] @ prev
                    KP = arc.left(s,prev)
                for t in range(d):
                    N[s,t] += phase*arc.join(arc.right(KX,t),y)
                    if not metric_only:
                        decorated = None
                        if contact:
                            decorated = contacts[s][t] if polynomial else transfer(A,apply_op(hh[s,:,:,t],A))
                        if separated:
                            inside = arc.right(KF,t) + KP @ Qright[t]
                            decorated = inside if decorated is None else inside+decorated
                        H[s,t] += phase*arc.join(decorated,y)
        del X, Y, y, F, prev, KX
        if not metric_only:
            del decorated, Pleft, Qright
            if L>2:
                del inside, KF, KP
            if polynomial:
                del contacts
        if metric_only:
            P = xp.linalg.matrix_power(E,L-1)
        else:
            prev,F = power_insert(E,V,L-2)
            P = prev @ E
        for s in range(d):
            N[s,s] += L*arc.diagonal(P)
        def assemble(block):
            return block.transpose(0,2,1,3).reshape(d*D*D,d*D*D)
        N = herm(assemble(N))
        if metric_only:
            return N
        # H currently contains one oriented arc only. Add half the coincident
        # term; the final Hermitian sum supplies the other half and other arc.
        for s in range(d):
            for t in range(d):
                left_edge = transfer(A, apply_op(hh[s,:,t,:], A))
                right_edge = transfer(A, apply_op(hh[:,s,:,t], A))
                diag = left_edge @ prev + prev @ right_edge
                if s == t:
                    diag += F
                H[s,t] += (L/2)*arc.diagonal(diag)
        H = assemble(H)
        return N, H + H.conj().T


def random_A(d, D, rng):
    A = rng.normal(size=(d,D,D)) + 1j*rng.normal(size=(d,D,D))
    return xp.asarray(A / np.sqrt(2*d*D))


def grow(A, D, rng):
    d, old, _ = A.shape
    if D < old:
        raise ValueError("bond schedule cannot shrink the loaded tensor")
    if D == old:
        return A
    out = .08*random_A(d,D,rng)
    out[:,:old,:old] = A
    return out


def balance(A):
    """Balance left/right Perron environments and return a local metric scale."""
    D = A.shape[1]
    def fixed(adj):
        e = xp.eye(D)/D
        for _ in range(256):
            f = sum(v.conj().T @ e @ v if adj else v @ e @ v.conj().T for v in A)
            f /= f.trace().real
            f = (f + e)/2  # also handles a near-period-two transfer
            if xp.linalg.norm(f-e) < 2e-13:
                break
            e = f
        return herm(f)
    def root(M):
        w, u = xp.linalg.eigh(M)
        return (u * xp.sqrt(xp.maximum(w,w[-1]*1e-12))) @ u.conj().T
    left, right = root(fixed(True)), root(fixed(False))
    _, lam, vh = xp.linalg.svd(left @ right)
    gauge = right @ (vh.conj().T / xp.sqrt(lam))
    A = xp.linalg.solve(xp.broadcast_to(gauge,A.shape), A @ gauge)
    weights = xp.sqrt(xp.maximum(lam, lam[0]*1e-5)/lam[0])
    return A, weights[:,None]*weights[None,:]


def optimize_ground(A, h, L, maxiter=350, rounds=3, ftol=2e-15):
    """L-BFGS in balanced gauges, finishing in a metric-preconditioned chart."""
    start = time.perf_counter()
    total = 0
    for iteration in range(rounds):
        A, S = balance(Ring(A,h,L).normalized())
        shape = A.shape
        Q = None
        if iteration>0 and iteration==rounds-1:
            N = Ring(A,h,L).tangent(0,metric_only=True)
            w,u = xp.linalg.eigh(N)
            keep = w>w[-1]*1e-12
            # Limit coordinate amplification during finite, nonlinear steps.
            Q = u[:,keep]/xp.sqrt(xp.maximum(w[keep],w[-1]*3e-9))
            del N, u
            initial = np.zeros(2*Q.shape[1])
        else:
            initial = np.ascontiguousarray(host(A*S)).view(float).ravel()

        def unpack(v):
            coeff = xp.asarray(v.view(np.complex128))
            if Q is None:
                return coeff.reshape(shape)/S
            return A+(Q@coeff).reshape(shape)

        def objective(v):
            a = unpack(v)
            # The derivative of this scalar cancels in the Rayleigh quotient.
            scale = float(xp.linalg.norm(a))/np.sqrt(shape[1])
            rr = Ring(a/scale,h,L)
            norm, e, nv, hv = rr.ground_env()
            g = (2*L/(norm*scale))*(hv-e*nv)
            g = g/S if Q is None else Q.conj().T@g.ravel()
            return e, np.ascontiguousarray(host(g)).view(float).ravel()

        result = opt.minimize(objective,initial,method="L-BFGS-B",jac=True,
                              options=dict(maxiter=maxiter,ftol=ftol,gtol=1e-11,
                                           maxls=45,maxcor=15))
        A = unpack(result.x)
        total += result.nit
        if result.nit < 3:
            break
    A, _ = balance(A)
    A = Ring(A,h,L).normalized()
    return A, dict(iterations=total, optimizer=str(result.message), seconds=time.perf_counter()-start)


def whitening(N, tol):
    ev,U = xp.linalg.eigh(N)
    e = host(ev)
    cut = float(tol*e[-1])
    k = int(np.searchsorted(e,max(0.,cut),side='right'))
    info = dict(relative_threshold=tol, threshold=cut, kept=len(e)-k,
                largest_discarded=float(e[k-1]) if k else None,
                smallest_kept=float(e[k]) if k<len(e) else None, min_metric=float(e[0]))
    return U[:,k:] / xp.sqrt(ev[k:]), info


def remove_reference(Q,N,ring):
    # A/L represents psi. Apply a complex Householder reflector implicitly.
    u = Q.conj().T @ (N @ (ring.A.ravel()/ring.L))
    length = float(xp.linalg.norm(u))
    if length==0:
        return Q
    u /= length
    u[0] += u[0]/abs(u[0]) if abs(u[0]) else 1
    u /= xp.linalg.norm(u)
    return Q[:,1:] - 2*(Q @ u)[:,None]*u[None,1:].conj()


def excitations(ring, m, E0, tol=1e-11, nev=12, cutoff=None):
    sync()
    start = time.perf_counter()
    N,H = ring.tangent(m,E0)
    sync()
    mid = time.perf_counter()
    Q,info = whitening(N,tol)
    if m % ring.L == 0:
        Q = remove_reference(Q,N,ring)
    del N
    heff = herm(Q.conj().T @ (H @ Q))
    del H
    ev, U = xp.linalg.eigh(heff)
    count = min(len(ev),nev) if cutoff is None else int(xp.sum(ev<=cutoff))
    B = Q @ U[:,:count]
    sync()
    info.update(build_seconds=mid-start, solve_seconds=time.perf_counter()-mid)
    return host(ev[:count]), B, info


def diagnostics(ring, tol):
    norm,e,nv,hv = ring.ground_env()
    N = ring.tangent(0,metric_only=True)
    Q,info = whitening(N,tol)
    grad = 2*ring.L*float(xp.linalg.norm(Q.conj().T @ (hv-e*nv).ravel()))/np.sqrt(norm)
    del N, Q
    return dict(E0=e, density=e/ring.L, norm=norm, gradient=grad,
                variance_per_site=ring.weight(0,e/ring.L)/ring.L, metric=info)


def local_dense(h, j, L, d):
    """Independent ED embedding, including the ordered bond (L-1,0)."""
    if not 0 <= j < L:
        raise ValueError("invalid site")
    if j == L-1:
        block = np.kron(h,np.eye(d**(L-2)))
        order = [*range(1,L),0]
        return block.reshape((d,)*(2*L)).transpose(order+[i+L for i in order]).reshape(block.shape)
    return np.kron(np.eye(d**j),np.kron(h,np.eye(d**(L-j-2))))


def explicit_state(A,L,B=None,m=0):
    """Literal trace product; deliberately independent of the transfer code."""
    d,D,_ = A.shape
    words = np.indices((d,)*L).reshape(L,-1)
    amps = np.zeros(d**L,complex)
    for j in range(1 if B is None else L):
        v = np.eye(D)
        for i,indices in enumerate(words):
            v = v @ (B if i==j and B is not None else A)[indices]
        amps += np.exp(-1j*TAU*m*j/L)*np.trace(v,axis1=-2,axis2=-1)
    return amps


def translate(v,L,d=2):
    return np.moveaxis(v.reshape((d,)*L),-1,0).ravel()


def brute_checks(verbose=True):
    worst = {}
    for L in (6,7):
        rng = np.random.default_rng(110+L)
        a,b,c = [host(random_A(2,3,rng)) for _ in range(3)]
        psi = explicit_state(a,L)
        norm = np.vdot(psi,psi).real
        def compare(label, actual, expected, scale=None):
            # Use vector norm for null identities; otherwise use relative error.
            den = np.linalg.norm(expected) if scale is None else scale
            err = float(np.linalg.norm(host(actual)-expected)/max(float(den),1e-300))
            worst[label] = max(worst.get(label,0.),err)
            if not err < 1e-10:
                raise AssertionError(f"L={L} {name} {label}: {err:.4g}, got {actual}, expected {expected}")
        for name in ('chir','xxz','tfim'):
            h = model(name)
            ring = Ring(a,h,L)
            terms = np.array([local_dense(host(h),j,L,2) for j in range(L)])
            H = terms.sum(axis=0)
            Hp = H @ psi
            ed = (np.vdot(psi,Hp)/norm).real
            norm_tm,e,nv,hv = ring.ground_env()
            compare('norm',norm_tm,norm)
            compare('energy numerator',e*norm_tm,ed*norm)
            compare('variance/site',ring.weight(0,e/L)/L,np.linalg.norm(Hp-ed*psi)**2/(L*norm))
            def energy(t,direction):
                v = explicit_state(a+t*direction,L)
                return (np.vdot(v,H@v)/np.vdot(v,v)).real
            def derivative(direction):
                step = 0.002
                # Eighth-order centered first derivative.
                return sum(q*(energy(k*step,direction)-energy(-k*step,direction))
                           for k,q in enumerate((0,4/5,-1/5,4/105,-1/280)) if k)/step
            gradient = 2*L*host(hv-e*nv)/norm_tm
            for direction in (b,1j*b,c):
                compare('gradient (directional)',np.vdot(direction,gradient).real,derivative(direction))
            derivatives = np.empty_like(a)
            for i in range(a.size):
                basis = np.zeros_like(a)
                basis.flat[i] = 1
                derivatives.flat[i] = derivative(basis)+1j*derivative(1j*basis)
            compare('gradient (all coordinates)',gradient,derivatives)
            for m in (0,1,-1,2,-2):
                p = TAU*m/L
                v,w = explicit_state(a,L,b,m),explicit_state(a,L,c,m)
                N,He = ring.tangent(m)
                N,He = host(N),host(He)
                if name=='chir':
                    Nk,Hk = ring.tangent(m,polynomial=True)
                    compare('polynomial N',Nk,N)
                    compare('polynomial H',Hk,He)
                    evals,cols,_ = excitations(ring,m,e,nev=4)
                    phys = projected_ritz(a,host(cols),L,m)
                    compare('Ritz norms',phys.conj().T@phys,np.eye(len(evals)))
                    if m==0:
                        compare('Ritz ground overlap',phys.conj().T@psi,0,np.linalg.norm(psi))
                compare('T psi',translate(psi,L),psi)
                compare('T Phi',translate(v,L),np.exp(1j*p)*v)
                compare('N',np.vdot(c,N@b.ravel()),np.vdot(w,v))
                compare('H',np.vdot(c,He@b.ravel()),np.vdot(w,H@v))
                X = rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
                # AX puts X after the defect; XA puts it one bond earlier.
                null = a@X - np.exp(1j*p)*(X@a)
                compare('gauge null',N@null.ravel(),0,np.linalg.norm(N)*np.linalg.norm(null))
                if m:
                    n = -m
                    ks = L/TAU*np.einsum('j,jab->ab',np.exp(1j*TAU*n*(np.arange(L)+.5)/L),terms)
                    target = ks @ psi/np.sqrt(norm)
                    compare('T Ht_n psi',translate(target,L),np.exp(1j*p)*target)
                    compare('M',np.vdot(b,host(ring.ks_vector(n))),np.vdot(v,target))
                    compare('W',ring.weight(n)*(L/TAU)**2,np.vdot(target,target))
    if verbose:
        print("Brute force, maximum errors",json.dumps(worst),flush=True)
    return worst


def momentum_basis(L,d,m):
    """Fourier transform each computational-basis translation orbit."""
    dim = d**L
    visited = np.zeros(dim,bool)
    orbits = []
    for initial in range(dim):
        if visited[initial]:
            continue
        orbit = []
        i = initial
        while not visited[i]:
            visited[i] = True
            orbit.append(i)
            i = (i%d)*(dim//d) + i//d
        if len(orbit)*m % L == 0:
            orbits.append(orbit)
    basis = np.zeros((dim,len(orbits)),complex)
    for j,orbit in enumerate(orbits):
        basis[orbit,j] = np.exp(-1j*TAU*m*np.arange(len(orbit))/L)/np.sqrt(len(orbit))
    return basis


def dense_ed(h,L,ms):
    """Full Hamiltonian from Kronecker products; solve in T eigenspaces."""
    if L > 12:
        raise ValueError("dense ED is limited to L <= 12")
    h = host(h)
    d = math.isqrt(len(h))
    H = local_dense(h,0,L,d)
    for j in range(1,L):
        H += local_dense(h,j,L,d)
    sectors = {}
    for m in ms:
        basis = momentum_basis(L,d,m)
        e,u = sla.eigh(herm(basis.conj().T @ (H @ basis)))
        sectors[m] = (e,basis@u)
    return H,sectors


def projected_ritz(A,B,L,m):
    """Explicit tangent vectors for overlap matching in the small ED check."""
    d,D,_ = A.shape
    words = np.indices((d,)*L).reshape(L,-1)
    v = np.eye(D)
    for site in words[1:]:
        v = v @ A[site]
    col = B.T.reshape(-1,d,D,D)
    chi = np.einsum('xvij,vji->vx',col[:,words[0]],v)
    result = np.zeros_like(chi)
    for j in range(L):
        result += np.exp(-1j*TAU*m*j/L)*chi
        chi = np.moveaxis(chi.reshape((d,)*L+(-1,)),L-1,0).reshape(chi.shape)
    return result


def groups(e,tol=1e-9):
    return np.split(np.arange(len(e)),np.flatnonzero(np.diff(e)>tol)+1)


def matched_weights(e,U,states,target):
    """Maximum-overlap assignment, with polar alignment inside ED multiplets."""
    overlap = U.conj().T @ states
    blocks = groups(e)
    scores = np.zeros(overlap.shape)
    for block in blocks:
        scores[block] = np.sum(abs(overlap[block])**2,axis=0)
    ranks, columns = opt.linear_sum_assignment(-scores)
    match = np.empty(states.shape[1],int)
    match[columns] = ranks
    weights = np.empty(len(match))
    fidelity = np.empty(len(match))
    for block in blocks:
        cols = np.flatnonzero(np.isin(match,block))
        if len(cols):
            proj = overlap[np.ix_(block,cols)]
            left,_,right = sla.svd(proj,full_matrices=False)
            aligned = U[:,block] @ (left@right)
            weights[cols] = abs(aligned.conj().T @ target)**2
            fidelity[cols] = np.abs(np.diag(aligned.conj().T @ states[:,cols]))**2
    return match,weights,fidelity


def ks_data(ring,n,ev,B):
    x = ev*ring.L/TAU
    M = host(B.conj().T @ ring.ks_vector(n))
    w = np.abs(M)**2
    W = ring.weight(n)*(ring.L/TAU)**2
    # At a lattice selection rule, W may vanish to roundoff.
    density = ring.h
    if n%ring.L:
        density = density-xp.trace(density)*xp.eye(len(density))/len(density)
    noise = 16*np.finfo(float).eps*(ring.L**2/TAU*float(xp.linalg.norm(density)))**2
    frac = float(w.sum()/W) if W>noise else None
    cs = [float(8*y/t**2) if t!=0 else None for y,t in zip(w,x)] if abs(n)==2 else None
    return dict(n=n,x=x.tolist(),M_real=M.real.tolist(),M_imag=M.imag.tolist(),
                weights=w.tolist(),W=W,weight_roundoff_estimate=noise,captured_fraction=frac,c=cs)


def refine_ground(A,h,L,tol=1e-11,steps=5):
    """Backtracked energy-Newton step in the orthogonal physical tangent.

    The curvature approximation is the centered tangent Hamiltonian. A small
    spectral floor handles near degeneracies and an indefinite early iterate.
    Using the one-hole residual avoids subtracting noisy whitened Ritz vectors.
    """
    count = 0
    for _ in range(steps):
        ring = Ring(A,h,L)
        _,e,nv,hv = ring.ground_env()
        N,H = ring.tangent(0,e)
        Q,_ = whitening(N,tol)
        Q = remove_reference(Q,N,ring)
        if Q.shape[1] == 0:
            break
        w,u = xp.linalg.eigh(herm(Q.conj().T @ (H @ Q)))
        f = Q.conj().T @ (L*(hv-e*nv)).ravel()
        floor = max(1e-7,0.05*float(xp.linalg.norm(f)))
        move = -Q @ (u @ ((u.conj().T @ f)/xp.maximum(w, floor)))
        del N, H, Q, u  # release dense matrices before backtracking or rebuilding
        for power in range(18):
            trial = A + (0.5**power)*move.reshape(A.shape)
            if L>=48:
                trial = Ring(trial,h,L).normalized()
            trial_ring = Ring(trial,h,L)
            candidate = trial_ring.ground_env()[1]
            if np.isfinite(candidate) and candidate < e - 1.5e-15*max(1,abs(e)):
                A = trial_ring.normalized()
                count += 1
                break
        else:
            break
    return A,count


def solve_ground(A,h,L,options):
    start = time.perf_counter()
    A,meta = optimize_ground(A,h,L,options['maxiter'],options['rounds'])
    A,n = refine_ground(A,h,L,options['tol'],options['refine'])
    ring = Ring(A,h,L)
    meta.update(diagnostics(ring,options['tol']))
    meta.update(D=A.shape[1],newton_steps=n,seconds=time.perf_counter()-start)
    return A,meta


def ed_checks():
    rows = []
    ms = (0,1,-1,2,-2)
    for name in ('chir','xxz','tfim'):
        L = 10
        h = model(name)
        _,sectors = dense_ed(h,L,ms)
        ebase = sectors[0][0][0]
        opts = dict(maxiter=350,rounds=3,refine=8,tol=1e-11)
        rng = np.random.default_rng(1)
        A = random_A(2,4,rng)
        print("ED",name,"reference",ebase,flush=True)
        for D in (4,8,12,16):
            A,meta = solve_ground(grow(A,D,rng),h,L,opts)
            ring = Ring(A,h,L)
            e0 = meta['E0']
            e,U = sectors[0]
            # Within an exactly degenerate p=0 ground space choose the
            # normalized projection of the optimized reference.
            Ug = U[:,groups(e)[0]]
            ref = Ug @ (Ug.conj().T @ explicit_state(host(A),L))
            ref /= np.linalg.norm(ref)
            row = dict(model=name,D=D,ed_reference=float(ebase),ground=meta,
                       ground_error=float(e0-ebase),sectors={})
            for m in ms:
                ev,B,info = excitations(ring,m,e0,nev=4)
                e,U = sectors[m]
                start = int(m==0)
                err = ev-(e[start:start+len(ev)]-ebase)
                sec = dict(omega=ev.tolist(),ed_omega=(e[start:start+len(ev)]-ebase).tolist(),
                           energy_errors=err.tolist(),common_reference_errors=(err+(e0-ebase)).tolist(),
                           metric=info)
                if m:
                    n = -m
                    target = np.zeros_like(ref)
                    for j in range(L):
                        target += L/TAU*np.exp(1j*TAU*n*(j+.5)/L)*(local_dense(host(h),j,L,2)@ref)
                    vecs = projected_ritz(host(A),host(B),L,m)
                    match,w,f = matched_weights(e,U,vecs,target)
                    data = ks_data(ring,n,ev,B)
                    sec.update(ks=data,ed_weights=w.tolist(),ed_ranks=(match+1).tolist(),
                               overlap_squared=f.tolist())
                row['sectors'][m] = sec
                print("ED",name,D,m,"dE0 %.5g"%(e0-ebase),"domega",np.array2string(err,precision=4),
                      "common",np.array2string(err+(e0-ebase),precision=4),
                      ("w "+np.array2string(np.array(sec['ks']['weights']),precision=4)
                       +" ED "+np.array2string(w,precision=4)+" ranks "+str(sec['ed_ranks'])) if m else "",
                      f"build={info['build_seconds']:.3f}s solve={info['solve_seconds']:.3f}s",
                      flush=True)
            print(f"Ground seconds={meta['seconds']:.3f} gradient={meta['gradient']:.4g} "
                  f"variance/L={meta['variance_per_site']:.4g}",flush=True)
            rows.append(row)
    return rows


def cft_checks():
    rows = []
    for name,Ds in (('tfim',(8,12)),('xxz',(8,12,16))):
        h = model(name,Jz=0.)
        for L in (16,32):
            rng = np.random.default_rng(1)
            A,_ = optimize_ground(random_A(2,4,rng),h,L)
            for D in Ds:
                opts = dict(maxiter=800,rounds=4,refine=6,tol=1e-11)
                A,meta = solve_ground(grow(A,D,rng),h,L,opts)
                ring = Ring(A,h,L)
                ms = (0,1,-1,2,-2) if name=="tfim" else (2,-2)
                row = dict(model=name,L=L,D=D,ground=meta,sectors={})
                for m in ms:
                    ev,B,info = excitations(ring,m,meta['E0'],nev=8)
                    sec = dict(x=(ev*L/TAU).tolist(),metric=info)
                    if abs(m)==2:
                        ks = ks_data(ring,-m,ev,B)
                        s = int(np.argmax(ks['weights']))
                        sec.update(ks=ks,stress_index=s,velocity=sec['x'][s]/2,
                                   stress_c=ks['c'][s])
                    row['sectors'][m] = sec
                v = np.mean([row['sectors'][m]['velocity'] for m in (2,-2)])
                row['velocity'] = float(v)
                for m,sec in row['sectors'].items():
                    sec['dimensions'] = (np.array(sec['x'])/v).tolist()
                print('CFT',name,L,D,'E0',meta['E0'],'var/L',meta['variance_per_site'],
                      'v',v,'c(+2,-2)',*[row['sectors'][m]['stress_c'] for m in (2,-2)],
                      'Delta', {m:np.round(s['dimensions'][:4],6).tolist()
                                for m,s in row['sectors'].items()},flush=True)
                print("ground seconds",meta['seconds'],"build/solve seconds",
                      {m:(s['metric']['build_seconds'],s['metric']['solve_seconds'])
                       for m,s in row['sectors'].items()},flush=True)
                rows.append(row)
    return rows


def timing_check():
    L,D = 64,24
    h = model('tfim')
    A = random_A(2,D,np.random.default_rng(64))
    ring = Ring(Ring(A,h,L).normalized(),h,L)
    e = ring.ground_env()[1]
    _,_,info = excitations(ring,2,e,nev=4)
    print("Timing L=64 D=24 m=2, normalized random tensor", info,flush=True)
    return info


def integer_list(value, L=None):
    if value == 'all' and L is not None:
        return list(range(L))
    if not value:
        return []
    values = []
    for item in value.split(','):
        if ':' in item:
            lo,hi = map(int,item.split(':'))
            step = 1 if hi >= lo else -1
            values.extend(range(lo,hi+step,step))
        else:
            values.append(int(item))
    return values


def gap_text(info):
    below = info['largest_discarded']
    return (f"tol={info['relative_threshold']:.2g} cut={info['threshold']:.4g} "
            f"kept={info['kept']} gap=({below if below is None else format(below,'.4g')}, "
            f"{info['smallest_kept'] if info['smallest_kept'] is None else format(info['smallest_kept'],'.4g')})")


def calculation(params):
    name,L = params['model'],params['L']
    h = model(name,params['delta'],params['h'],params['Jz'],params['file'])
    d = math.isqrt(len(h))
    ds = params['Dsched']
    rng = np.random.default_rng(params['seed'])
    if params['load']:
        A = xp.asarray(np.load(params['load'],allow_pickle=False),dtype=xp.complex128)
        if A.ndim!=3 or A.shape[0]!=d or A.shape[1]!=A.shape[2] or not xp.all(xp.isfinite(A)):
            raise ValueError("loaded A must be finite with shape (d,D,D)")
        peak = float(xp.max(xp.abs(A)))
        if not peak:
            raise ValueError("zero tensor in load file")
        A = Ring(A/peak,h,L).normalized()
        if ds is None:
            ds = [A.shape[1]]
    else:
        ds = ds or [4,8]
        A = random_A(d,ds[0],rng)
    if not ds or min(ds)<1 or any(y<x for x,y in zip(ds,ds[1:])):
        raise ValueError("Dsched must be a nondecreasing list of positive dimensions")
    params['Dsched'] = ds
    ms = {}
    for m in params['ps']+[-n for n in params['ns']]:
        ms.setdefault(m%L,m)
    params['computed_momenta'] = list(ms.values())
    print("Parameters",json.dumps(params),flush=True)
    out = dict(parameters=params,ground=[],spectra=[],koo_saleur=[])
    for D in ds:
        A,ground = solve_ground(grow(A,D,rng),h,L,params)
        out['ground'].append(ground)
        print(f"D={D} E0={ground['E0']:.15g} E0/L={ground['density']:.13g} "
              f"norm={ground['norm']:.12g} |grad E|={ground['gradient']:.5g} "
              f"var/L={ground['variance_per_site']:.6g} seconds={ground['seconds']:.3f} "
              f"iterations={ground['iterations']} Newton={ground['newton_steps']} "
              f"L-BFGS={ground['optimizer']}",flush=True)
    e = ground['E0']
    ring = Ring(A,h,L)
    ks = {}
    for m in ms.values():
        ev,B,info = excitations(ring,m,e,params['tol'],params['nev'],params['cutoff'])
        sec = dict(m=m,omega=ev.tolist(),x=(ev*L/TAU).tolist(),negative=(ev<0).tolist(),**info)
        out['spectra'].append(sec)
        print(f"m={m} {gap_text(info)} build={info['build_seconds']:.3f}s "
              f"solve={info['solve_seconds']:.3f}s\nm s omega x",flush=True)
        for i,(v,x) in enumerate(zip(ev,sec['x']),1):
            print(f"{m:4d} {i:4d} {v:17.11g} {x:17.11g}"+(" NEGATIVE" if v<0 else ""))
        for n in params['ns']:
            if (m+n)%L==0:
                ks[n] = ks_data(ring,n,ev,B)
                ks[n]['m'] = m
        del B
    for n in params['ns']:
        sec = ks[n]
        out['koo_saleur'].append(sec)
        frac = sec['captured_fraction']
        print(f"n={n} m={sec['m']} W={sec['W']:.12g} captured_fraction="
              + (str(frac) if frac is not None else 'undefined (weight at roundoff)'))
        print("s x Re(M) Im(M) |M|^2"+(" c" if sec['c'] is not None else ""))
        for i,(x,re,im,w) in enumerate(zip(sec['x'],sec['M_real'],sec['M_imag'],sec['weights'])):
            c = sec['c']
            print(f"{i+1:4d} {x:15.9g} {re:15.9g} {im:15.9g} {w:15.9g}"
                  +(f" {c[i]}" if c is not None else ""))
    return A,out


def self_checks():
    start = time.perf_counter()
    brute = brute_checks()
    ed = ed_checks()
    cft = cft_checks()
    timing = timing_check()
    for name in ('chir','xxz','tfim'):
        subset = [r for r in ed if r['model']==name]
        def errors(row):
            energy = max(abs(t) for s in row['sectors'].values() for t in s['energy_errors'])
            weight = max(abs(a-b) for s in row['sectors'].values() if 'ks' in s
                         for a,b in zip(s['ks']['weights'],s['ed_weights']))
            return np.array([abs(row['ground_error']),energy,weight])
        first,last = errors(subset[0]), errors(subset[-1])
        if not (np.all(last < first) and np.all(last < [1e-5,1e-4,2e-3])):
            raise AssertionError(f"ED {name}: initial/final errors {first}, {last}")
    for row in cft:
        if row['model']=='tfim':
            dims = row['sectors'][0]['dimensions']
            if max(abs(dims[0]-.125),abs(dims[1]-1))>.03:
                raise AssertionError("TFIM scaling dimensions outside finite-size check")
        target = .5 if row['model']=='tfim' else 1.
        for m in (-2,2):
            if row['D']==(12 if target==.5 else 16):
                if abs(row['sectors'][m]['stress_c']-target)>.12:
                    raise AssertionError("large-D stress weight outside finite-size check")
    elapsed = time.perf_counter()-start
    print(f"Self-checks passed in {elapsed:.2f}s",flush=True)
    return dict(brute=brute,ed=ed,cft=cft,timing=timing,seconds=elapsed)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Options are key=value. Example:
  python3 pumps.py model=chir L=64 Dsched=8,16,24,32 ps=-4:4 ns=1,2 nev=30 save=out/chir
  python3 pumps.py check=1 save=out/checks

model=chir|xxz|tfim|file  delta=.1875 h=1.783 Jz=.5 file=path.npy
L=16 Dsched=4,8 load=path.npy seed=1 ps=-2:2|all ns=1,2
nev=12 cutoff=energy (overrides nev) tol=1e-11
maxiter=350 (per L-BFGS round) rounds=3 refine=5 (tangent Newton steps)
save=pumps_out (prefix for .json and _A.npy)
check=1 alone runs the numerical suite. Add model parameters to also run a fit.
PUMPS_BACKEND=numpy|cupy selects the array library before import.""")
    parser.add_argument('options',nargs='*',metavar='key=value')
    items = parser.parse_args(argv).options
    defaults = dict(model='chir',L=16,Dsched=None,ps='-2:2',ns='1,2',nev=12,
                    delta=.1875,h=1.783,Jz=.5,file=None,seed=1,save='pumps_out',load=None,
                    check=0,cutoff=None,tol=1e-11,maxiter=350,rounds=3,refine=5)
    options = {}
    for item in items:
        key,sep,value = item.partition('=')
        if not sep or key not in defaults:
            parser.error(f"unknown argument {item}")
        options[key] = value
    params = {**defaults,**options}
    try:
        for k in ('L','nev','seed','check','maxiter','rounds','refine'):
            params[k] = int(params[k])
        for k in ('delta','h','Jz','tol','cutoff'):
            if params[k] is not None:
                params[k] = float(params[k])
        L = params['L']
        if (L<2 or params['nev']<0 or not 0<params['tol']<1 or params['check'] not in (0,1)
            or params['maxiter']<1 or params['rounds']<1 or params['refine']<0):
            raise ValueError("invalid L, nev, tol, check or optimizer limit")
        if params['Dsched'] is not None:
            params['Dsched'] = integer_list(params['Dsched'])
            if not params['Dsched'] or min(params['Dsched'])<1:
                raise ValueError("Dsched requires positive integers")
        params['ps'] = integer_list(params['ps'],L)
        params['ns'] = list(dict.fromkeys(integer_list(params['ns'])))
    except ValueError as exc:
        parser.error(str(exc))
    params['backend'] = xp.__name__
    prefix = Path(params['save']).expanduser()
    prefix.parent.mkdir(parents=True,exist_ok=True)
    report = {}
    if not params['check'] or options.keys()-{'check','save'}:
        A,report = calculation(params)
        np.save(str(prefix)+'_A.npy',host(A),allow_pickle=False)
    if params['check']:
        print("Self-check parameters",json.dumps(dict(backend=xp.__name__,check=1,save=str(prefix))),
              flush=True)
        report['self_checks'] = self_checks()
    path = Path(str(prefix)+'.json')
    path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print("Saved",path,flush=True)


if __name__=='__main__':
    main()
