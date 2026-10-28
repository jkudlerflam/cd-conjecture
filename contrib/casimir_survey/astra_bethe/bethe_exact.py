#!/usr/bin/env python3
"""Full finite-size nested Bethe equations, NumPy and SciPy only.

Usage
=====
python3 bethe_exact.py --validate-only
python3 bethe_exact.py
python3 bethe_exact.py --families 2 3 --max-multiple 60

The default run validates EVERY supplied ED energy before production, then
continues k=2 solutions in L in each requested family. Each family has a
20-minute wall-clock budget, checked between Newton iterations. Default size
ceilings are 30N for N<=5 and 15N for N=7,8; override --max-multiple, or use
--extend to continue until the budget is exhausted. --resume reuses verified
production roots from the output JSON (validation is still rerun).

Equations and Newton derivative
==============================
e_n(z)=(z+in/2)/(z-in/2). With M_r=k L(N-r)/N, the original equations are
 e_k(z_ra)^(L delta_r1) = prod_{b != a} e_2(z_ra-z_rb)
                         prod_{s=r+-1,b} e_-1(z_ra-z_sb).
Thus the principal logarithmic residual is the imaginary-part-wrapped sum
 F_ra=L delta_r1 log e_k(z_ra) - sum_{b != a} log e_2(z_ra-z_rb)
                                      + sum_{s=r+-1,b} log e_1(z_ra-z_sb).
Writing D_n(z)=1/(z+in/2)-1/(z-in/2), its diagonal derivative is
 L delta_r1 D_k - sum_{b != a} D_2 + sum_{s=r+-1,b} D_1.
Off-diagonal same-level entries are +D_2, adjacent-level entries are -D_1.
All other entries vanish. Locally the branch wrapping has zero derivative.
Newton solves J dz=-F, with backtracking in the norm of the full residual.

Strings enter ONLY as seeds. No fused or ideal-string equation is used for
the final roots or energy. The direct log formulation suffices for the tested
ground states: their string deviations remain resolved in double precision.
The script does not claim to resolve exponentially small deviations, singular
roots, repeated roots, or roots at infinity. Such configurations are rejected,
and a family stops if all seeds fail. For k=2, conjugation and reflection are
enforced in a real linear coordinate subspace. Both x and y in x+-iy vary.

Independent checks and energy
=============================
An independent evaluation of the ORIGINAL products checks
each root's abs(log(LHS/RHS)), without any fusion or cancellation of internal
string scattering. It uses NumPy clongdouble, which has extra precision only
on platforms that support it (on Apple arm64 it is the same as complex128).
E=L H_k - (1/2) sum_{level 1} k/(z*z+k*k/4).
We check its imaginary part, conjugation, reflection, root count, distinctness,
and finite root size. All roots, each member residual, and rejected attempts
are saved. Small residuals alone do not imply a physical Bethe state.
ED selects no fitted parameter: candidate acceptance uses only equations and
root diagnostics. An ED mismatch then fails validation and blocks production.
All validated ground states retain deformed, self-conjugate k-string groupings.
At larger sizes the ground-state assignment relies on continuation from those
states; for N=5,7,8 no independent ED targets were supplied.

Fits
====
The analytic bulk energy is used only AFTER solving the roots. Fits use
 A(L)=-6L(E-e_bulk L)/pi = A_inf + a/L^2 [+ b/L^4],
or include 1/log(L), 1/log(L)^3, and optionally 1/log(L)^4 corrections.
Free-bulk fits instead fit E/L=e_fit-pi/(6L^2)*(A_inf+corrections).
Columns are normalized before linear least squares. Several late-size windows
are used, with at least four degrees of freedom. The reported central value is
the median over fits and the uncertainty is the maximum displacement from it.
This is fit-form/window sensitivity, not a statistical confidence interval or
a controlled bound on omitted logarithms. Every individual fit is retained.
"""
import argparse
import json
import os
import time
import warnings
from pathlib import Path

# Small/medium dense systems usually benefit from avoiding thread oversubscription.
# Honor thread settings explicitly supplied by the caller.
for _name in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):
    os.environ.setdefault(_name,'1')
import numpy as np
from scipy.linalg import solve, LinAlgWarning
from scipy.special import digamma

DEADLINE = float('inf')

def check_deadline():
    if time.monotonic() > DEADLINE:
        raise TimeoutError('Per-family time budget exhausted')

TARGETS = {
    (2, 2): [(4, 1.6492189406), (8, 3.8393150668), (12, 5.8948058277)],
    (2, 3): [(8, 4.9267655562), (12, 7.5551867405)],
    (2, 4): [(8, 5.7817491256), (10, 7.3280616238)],
    (3, 2): [(6, 1.2020467253), (9, 2.0739788158)],
    (3, 3): [(6, 1.5976165898)],
    (4, 2): [(8, 0.9558401418)],
}

def phase(x, terms, derivative=False):
    out = np.zeros_like(x, dtype=float)
    for n, weight in terms:
        out += weight*(n/(x*x+n*n/4) if derivative else 2*np.arctan(2*x/n))
    return out

def initial_centres(N, L):
    blocks = []
    for r in range(1, N):
        n = L*(N-r)//N
        u = (np.arange(n)+0.5)/L-(N-r)/(2*N)
        blocks.append(N/np.pi*np.arctanh(np.tan(np.pi*u)*np.tan(np.pi*r/(2*N))))
    return np.concatenate(blocks)

def slices_for(N, k, L):
    ends = np.cumsum([0]+[k*L*(N-r)//N for r in range(1, N)])
    return [slice(ends[r], ends[r+1]) for r in range(N-1)]

def solve_centres(N, k, L):
    slices = slices_for(N, 1, L)
    source = [(n, 1) for n in range(1, 2*k, 2)]
    same = [(n, 2) for n in range(2, 2*k, 2)]+[(2*k, 1)]
    quantum = np.concatenate([np.arange(s.stop-s.start)-(s.stop-s.start-1)/2 for s in slices])
    def eq(x):
        f = -2*np.pi*quantum.copy()
        J = np.zeros((len(x), len(x)))
        for r, sl in enumerate(slices):
            xr = x[sl]
            diagonal = np.zeros(len(xr))
            if r == 0:
                f[sl] += L*phase(xr, source)
                diagonal += L*phase(xr, source, True)
            d = xr[:, None]-xr[None, :]
            f[sl] -= phase(d, same).sum(axis=1)
            kernel = phase(d, same, True)
            np.fill_diagonal(kernel, 0)
            diagonal -= kernel.sum(axis=1)
            J[sl, sl] = kernel
            for s in (r-1, r+1):
                if 0 <= s < N-1:
                    ss = slices[s]
                    d = xr[:, None]-x[ss][None, :]
                    f[sl] += phase(d, source).sum(axis=1)
                    kernel = phase(d, source, True)
                    diagonal += kernel.sum(axis=1)
                    J[sl, ss] = -kernel
            indices = np.arange(sl.start, sl.stop)
            J[indices, indices] += diagonal
        return f, J
    x = initial_centres(N, L)
    for _ in range(80):
        check_deadline()
        f, J = eq(x)
        if np.max(np.abs(f)) < 1e-12:
            return x
        step = solve(J, -f, assume_a='sym')
        for i in range(40):
            trial = x+2.**(-i)*step
            if all(np.all(np.diff(trial[s])>0) for s in slices) and np.linalg.norm(eq(trial)[0]) < np.linalg.norm(f):
                x = trial
                break
        else:
            raise RuntimeError('Ideal-centre Newton stalled')
    raise RuntimeError('Ideal-centre Newton iteration limit')

def loge(n, z):
    return np.log(z+0.5j*n)-np.log(z-0.5j*n)

def dloge(n, z):
    return 1/(z+0.5j*n)-1/(z-0.5j*n)

def equations(z, N, k, L, jac=True):
    slices = slices_for(N, k, L)
    f = np.zeros(len(z), complex)
    J = np.zeros((len(z), len(z)), complex) if jac else None
    for r, sl in enumerate(slices):
        zr = z[sl]
        diagonal = np.zeros(len(zr), complex)
        if r == 0:
            f[sl] += L*loge(k, zr)
            if jac:
                diagonal += L*dloge(k, zr)
        d = zr[:, None]-zr[None, :]
        kernel = loge(2, d)
        np.fill_diagonal(kernel, 0)
        f[sl] -= kernel.sum(axis=1)
        if jac:
            kernel = dloge(2, d)
            np.fill_diagonal(kernel, 0)
            diagonal -= kernel.sum(axis=1)
            J[sl, sl] = kernel
        for s in (r-1, r+1):
            if 0 <= s < N-1:
                ss = slices[s]
                d = zr[:, None]-z[ss][None, :]
                f[sl] += loge(1, d).sum(axis=1)
                if jac:
                    kernel = dloge(1, d)
                    diagonal += kernel.sum(axis=1)
                    J[sl, ss] = -kernel
        if jac:
            indices = np.arange(sl.start, sl.stop)
            J[indices, indices] += diagonal
    f.imag = (f.imag+np.pi)%(2*np.pi)-np.pi
    return f, J

def newton(z, N, k, L, tol=2e-12, maxiter=100):
    z = np.array(z, complex)
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('error', LinAlgWarning)
        for it in range(maxiter):
            check_deadline()
            f, J = equations(z, N, k, L)
            residual = np.max(np.abs(f))
            if not np.isfinite(residual):
                return z, float('inf'), it
            if residual < tol:
                return z, residual, it
            try:
                step = solve(J, -f)
            except (np.linalg.LinAlgError, LinAlgWarning, ValueError):
                return z, residual, it
            for i in range(24):
                trial = z+2.**(-i)*step
                ff = equations(trial, N, k, L, False)[0]
                if np.linalg.norm(ff) < np.linalg.norm(f):
                    z = trial
                    break
            else:
                return z, residual, it
    return z, float(np.max(np.abs(equations(z, N, k, L, False)[0]))), maxiter

def diagnostics(z, N, k, L):
    slices = slices_for(N, k, L)
    zz = np.asarray(z[slices[0]],dtype=np.clongdouble)
    Hk = np.sum(1/np.arange(1,k+1,dtype=np.longdouble))
    energy = L*Hk-np.longdouble(.5)*np.sum(k/(zz**2+np.longdouble(k*k)/4))
    conj = refl = 0.
    sep = float('inf')
    for sl in slices:
        zr = z[sl]
        conj = max(conj, np.max(np.min(abs(zr[:, None]-zr.conj()[None, :]), axis=1)))
        refl = max(refl, np.max(np.min(abs(zr[:, None]+zr[None, :]), axis=1)))
        d = abs(zr[:, None]-zr[None, :])
        np.fill_diagonal(d, np.inf)
        sep = min(sep, np.min(d))
    return dict(E=float(energy.real), energy_imag=float(energy.imag),
                conjugation_error=float(conj), reflection_error=float(refl),
                minimum_separation=float(sep), maximum_root=float(np.max(abs(z))))

def product_check(z, N, k, L):
    """Independent check of ORIGINAL products, using clongdouble where available.

    No fused identities or string equations enter here. This is deliberately
    separate from the Newton residual. Report the principal log(LHS/RHS).
    """
    z = np.asarray(z, dtype=np.clongdouble)
    blocks = [z[s] for s in slices_for(N, k, L)]
    out = []
    def e(n, w):
        return (w+np.clongdouble(0.5j*n))/(w-np.clongdouble(0.5j*n))
    with np.errstate(all='ignore'):
        for r, zr in enumerate(blocks):
            vals = []
            for a, w in enumerate(zr):
                lhs = e(k, w)**L if r == 0 else np.clongdouble(1)
                factors = [e(2, w-np.delete(zr, a))]
                for s in (r-1, r+1):
                    if 0 <= s < N-1:
                        factors.append(e(-1, w-blocks[s]))
                factors = np.concatenate(factors)
                rhs = np.prod(factors)
                value = np.log(lhs/rhs)
                if not np.isfinite(value):
                    # Same original products, separated into magnitude and
                    # unit phase to prevent overflow at very large L.
                    src = e(k,w) if r==0 else np.clongdouble(1)
                    power = L if r==0 else 0
                    logmag = power*np.log(abs(src))-np.sum(np.log(abs(factors)))
                    unit = (src/abs(src))**power/np.prod(factors/abs(factors))
                    value = logmag+1j*np.angle(unit)
                vals.append(float(abs(value)))
            out.append(vals)
    return out

def symmetric_map(N, L):
    """k=2 conjugation/reflection fixed subspace. No fixed string spacing.

    Left-half roots are x_a +- i y_a, right-half roots -x_a +- i y_a.
    An odd number of pairs adds roots +-i y_0. Both x and y are free.
    B maps real independent coordinates to all complex roots. Selecting one
    member per symmetry orbit gives the same number of real equations.
    """
    nvar = L*(N-1)//2
    nroot = 2*nvar
    B = np.zeros((nroot, nvar), complex)
    picks, parts = [], []
    offset = column = 0
    for r in range(1, N):
        n = L*(N-r)//N
        m = n//2
        for a in range(m):
            b = n-1-a
            B[offset+2*a:offset+2*a+2, column+a] = 1
            B[offset+2*b:offset+2*b+2, column+a] = -1
            B[offset+2*a:offset+2*a+2, column+m+a] = [1j, -1j]
            B[offset+2*b:offset+2*b+2, column+m+a] = [1j, -1j]
            picks.extend([offset+2*a, offset+2*a])
            parts.extend([0, 1])
        if n%2:
            B[offset+2*m:offset+2*m+2, column+2*m] = [1j, -1j]
            picks.append(offset+2*m)
            parts.append(0)
        offset += 2*n
        column += n
    assert column == nvar and offset == nroot
    return B, np.array(picks), np.array(parts)

def symmetric_newton(z, N, L, tol=2e-12, maxiter=60):
    B, picks, parts = symmetric_map(N, L)
    # Orthogonal columns make projection explicit.
    t = (B.conj().T@z).real/np.sum(abs(B)**2, axis=0)
    columns = [(np.flatnonzero(B[:, j]), B[B[:, j]!=0, j]) for j in range(B.shape[1])]
    def select(f):
        ff = f[picks]
        return np.where(parts == 0, ff.real, ff.imag)
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('error', LinAlgWarning)
        for it in range(maxiter):
            check_deadline()
            z = B@t
            f, J = equations(z, N, 2, L)
            residual = np.max(abs(f))
            if not np.isfinite(residual) or residual < tol:
                return z, float(residual), it
            # Every column has at most four nonzero entries.
            selected = J[picks]
            jj = np.column_stack([selected[:, idx]@weights for idx, weights in columns])
            jr = np.where(parts[:, None] == 0, jj.real, jj.imag)
            ff = select(f)
            try:
                step = solve(jr, -ff)
            except (np.linalg.LinAlgError, LinAlgWarning, ValueError):
                return z, float(residual), it
            for i in range(22):
                trial = B@(t+2.**(-i)*step)
                ft = equations(trial, N, 2, L, False)[0]
                if np.linalg.norm(ft) < np.linalg.norm(f):
                    t += 2.**(-i)*step
                    break
            else:
                return z, float(residual), it
    return B@t, float(np.max(abs(equations(B@t,N,2,L,False)[0]))), maxiter

def continuation_seed(N, k, L, old_L, old_z):
    """Interpolate corrections on the string quantile grid, respecting parity.

    The ideal centres are solved afresh. Corrections are rescaled by old_L/L,
    a starting guess only. No interpolation constraint enters the equations.
    """
    old_x = solve_centres(N, k, old_L)
    new_x = solve_centres(N, k, L)
    old_slices = slices_for(N, 1, old_L)
    new_slices = slices_for(N, 1, L)
    shifts = .5j*(k+1-2*np.arange(1,k+1))
    delta = old_z.reshape(-1,k)-old_x[:,None]-shifts
    new = new_x[:,None]+np.broadcast_to(shifts,(len(new_x),k)).copy()
    for so, sn in zip(old_slices,new_slices):
        no, nn = so.stop-so.start, sn.stop-sn.start
        uo, un = (np.arange(no)+.5)/no, (np.arange(nn)+.5)/nn
        for j in range(k):
            new[sn,j] += old_L/L*np.interp(un,uo,delta[so,j])
    return new.ravel()

def seed(N, k, L, scale=0.05):
    x = solve_centres(N, k, L)
    shifts = 0.5j*(k+1-2*np.arange(1, k+1))
    return (x[:, None]+shifts[None, :]*(1+scale)).ravel()

def solve_ground(N,k,L,previous=None):
    if N<2 or k<1 or L<N or L%N:
        raise ValueError('Require N>=2, k>=1, L>=N, and L divisible by N')
    started = time.monotonic()
    guesses = []
    if previous is not None:
        old_L,old_z = previous
        guesses.append(('continuation',continuation_seed(N,k,L,old_L,old_z)))
    for scale in (.1,.03,.2,-.05,.5,.01):
        guesses.append((f'ideal_seed_scale={scale}',None))
    attempts = []
    for name,initial in guesses:
        if initial is None:
            initial = seed(N,k,L,float(name.split('=')[1]))
        z,residual,iterations = (symmetric_newton(initial,N,L) if k==2
                                 else newton(initial,N,k,L))
        d = diagnostics(z,N,k,L)
        attempt = dict(seed=name,newton_residual=float(residual),iterations=iterations,**d)
        attempts.append(attempt)
        valid = (np.isfinite(residual) and residual<2e-9
                 and all(np.isfinite(v) for v in d.values())
                 and abs(d['energy_imag'])<1e-9
                 and d['conjugation_error']<1e-8 and d['reflection_error']<1e-8
                 and d['minimum_separation']>1e-7
                 and d['maximum_root']<max(20.,10*N*np.log(L+1)))
        if valid:
            residuals = product_check(z,N,k,L)
            original = max(map(max,residuals))
            attempt['product_residual'] = original
            valid = np.isfinite(original) and original<2e-9
        attempt['accepted'] = bool(valid)
        if valid:
            pairs=z.reshape(-1,k)
            shifts=.5j*(k+1-2*np.arange(1,k+1))
            deviations=pairs-pairs.real.mean(axis=1)[:,None]-shifts
            return z,dict(N=N,k=k,L=L,**d,
                max_log_product_residual=original,log_product_residuals_by_level=residuals,
                roots_by_level=[[[float(v.real),float(v.imag)] for v in z[s]]
                                for s in slices_for(N,k,L)],
                max_string_deviation=float(np.max(abs(deviations))),
                elapsed_seconds=time.monotonic()-started,attempts=attempts)
    raise RuntimeError(f'No acceptable full-root solution at {(N,k,L)}. Attempts: {attempts}')

def bulk_energy(N,k):
    j=np.arange(1,k+1,dtype=float)
    return float(np.sum(1/j)-np.sum(digamma(1+(j-1)/N)-digamma(j/N))/N)

def fit_family(rows,N,k=2):
    if len(rows)<8:
        return dict(status='insufficient_sizes',fits=[])
    L=np.array([r['L'] for r in rows],float)
    E=np.array([r['E'] for r in rows],float)
    e0=bulk_energy(N,k)
    A=-6*L*(E-e0*L)/np.pi
    models={
        'L2': ['L2'], 'L2_L4': ['L2','L4'],
        'log1_L2': ['log1','L2'], 'log3_L2': ['log3','L2'],
        'log3_log4_L2': ['log3','log4','L2'],
    }
    starts=sorted(set([4*N,6*N,8*N]+[int(f*L[-1]/N)*N for f in (.35,.5,.65)]))
    fits=[]
    for start in starts:
        mask=L>=start
        ll=L[mask]
        for model,terms in models.items():
            cols=[np.ones(len(ll))]
            for term in terms:
                cols.append(ll**(-int(term[1:])) if term.startswith('L')
                            else np.log(ll)**(-int(term[3:])))
            X=np.column_stack(cols)
            for free in (False,True):
                if free:
                    xx=np.column_stack([np.ones(len(ll)),-np.pi/(6*ll[:,None]**2)*X])
                    yy=E[mask]/ll
                    labels=['e0','A_inf']+terms
                else:
                    xx=X
                    yy=A[mask]
                    labels=['A_inf']+terms
                if len(ll)<xx.shape[1]+4:
                    continue
                scales=np.linalg.norm(xx,axis=0)
                xn=xx/scales
                coef,_,rank,sing=np.linalg.lstsq(xn,yy,rcond=None)
                if rank<xx.shape[1] or sing[0]/sing[-1]>1e10:
                    continue
                coef/=scales
                coeff={label:float(v) for label,v in zip(labels,coef)}
                fits.append(dict(model=model,free_bulk=free,L_min=int(ll[0]),
                    L_max=int(ll[-1]),n_sizes=len(ll),coefficients=coeff,
                    rms=float(np.sqrt(np.mean((xx@coef-yy)**2))),
                    scaled_condition=float(sing[0]/sing[-1])))
    def spread(values):
        if not values:
            return None
        center=float(np.median(values))
        return dict(value=center,form_window_error=float(np.max(abs(np.array(values)-center))),
                    minimum=float(min(values)),maximum=float(max(values)),count=len(values))
    fixed=spread([f['coefficients']['A_inf'] for f in fits if not f['free_bulk']])
    freeA=spread([f['coefficients']['A_inf'] for f in fits if f['free_bulk']])
    freeE=spread([f['coefficients']['e0'] for f in fits if f['free_bulk']])
    Hk=float(np.sum(1/np.arange(1,k+1)))
    prediction=np.pi*k*(N*N-1)/(N*(N+k))
    return dict(status='fitted',bulk_formula=e0,A_prediction=prediction,
        norm_h=Hk/2,ratio_prediction=prediction/(np.pi*Hk),
        A_inf=fixed,A_inf_free_bulk=freeA,e0_free_fit=freeE,
        e0_free_minus_formula=None if freeE is None else freeE['value']-e0,
        normalized_A=None if fixed is None else dict(value=fixed['value']/(np.pi*Hk),
                        form_window_error=fixed['form_window_error']/(np.pi*Hk)),
        uncertainty_definition='Median and maximal displacement across fit forms and windows; not statistical.',
        fits=fits)

def write_results(path,result):
    """Atomic checkpoint; nonfinite failed-attempt diagnostics become null."""
    def clean(value):
        if isinstance(value,dict):
            return {str(k):clean(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):
            return [clean(v) for v in value]
        if isinstance(value,(float,np.floating)) and not np.isfinite(value):
            return None
        return value
    path=Path(path)
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(clean(result),separators=(',',':'),allow_nan=False)+'\n')
    temporary.replace(path)

def validation():
    rows=[]
    print('Validation against ED, using full complex roots',flush=True)
    print(' N k   L           E_Bethe            E_ED        delta       max |log ratio|',flush=True)
    for (N,k),targets in TARGETS.items():
        for L,target in targets:
            try:
                z,row=solve_ground(N,k,L)
                row.update(ED=target,delta=row['E']-target,passed=abs(row['E']-target)<1e-8)
                print(f'{N:2} {k:1} {L:3} {row["E"]:18.12f} {target:15.10f}'
                      f' {row["delta"]:+12.3e} {row["max_log_product_residual"]:12.3e}',flush=True)
            except RuntimeError as exc:
                row=dict(N=N,k=k,L=L,ED=target,passed=False,error=str(exc))
                print(str(exc),flush=True)
            rows.append(row)
    return rows

def print_summary(result):
    print('\n N k Lmax     A_inf (form/window spread)       prediction      A/(2 pi ||h||)        ratio prediction',flush=True)
    for family in result.get('families',[]):
        fit=family['fit']
        if fit.get('A_inf') is None:
            print(f'{family["N"]:2} 2   insufficient sizes for fits',flush=True)
            continue
        a,ratio=fit['A_inf'],fit['normalized_A']
        print(f'{family["N"]:2} 2 {family["data"][-1]["L"]:4}'
              f'  {a["value"]:12.8f} +/- {a["form_window_error"]:.3g}'
              f'     {fit["A_prediction"]:12.8f}'
              f'    {ratio["value"]:.7f} +/- {ratio["form_window_error"]:.2g}'
              f'    {fit["ratio_prediction"]:.7f}',flush=True)
    print('\n N      e0 formula          e0 free fit       form/window spread   difference',flush=True)
    for family in result.get('families',[]):
        fit=family['fit']
        free=fit.get('e0_free_fit')
        if free is not None:
            print(f'{family["N"]:2} {fit["bulk_formula"]:19.12f} {free["value"]:19.12f}'
                  f'   +/- {free["form_window_error"]:.3e}'
                  f'     {fit["e0_free_minus_formula"]:+.3e}',flush=True)
    print('Errors describe fit-form/window sensitivity, not statistical confidence intervals.',flush=True)

def main():
    global DEADLINE
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--validate-only',action='store_true')
    parser.add_argument('--families',type=int,nargs='+',default=[2,3,4,5,7,8])
    parser.add_argument('--max-multiple',type=int)
    parser.add_argument('--extend',action='store_true',help='Continue beyond default size ceilings until time budget')
    parser.add_argument('--minutes-per-family',type=float,default=20.)
    parser.add_argument('--output',default='bethe_exact_results.json')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if any(N<2 for N in args.families) or args.minutes_per_family<=0:
        parser.error('Require N>=2 and a positive time budget')
    if args.max_multiple is not None and args.max_multiple<2:
        parser.error('--max-multiple must be at least 2')
    saved=json.loads(Path(args.output).read_text()) if args.resume and Path(args.output).exists() else {}
    retained=[f for f in saved.get('families',[]) if f['N'] not in args.families]
    result=dict(schema_version=1,status='validating',validation=validation(),families=retained,
                configuration=vars(args),numpy_version=np.__version__,
                longdouble_epsilon=float(np.finfo(np.longdouble).eps))
    passed=all(row['passed'] for row in result['validation'])
    result['status']='validation_passed' if passed else 'validation_failed_no_production'
    write_results(args.output,result)
    if not passed:
        print('VALIDATION FAILED. Production and extrapolation were not run.',flush=True)
        return 1
    print('All 11 ED validations passed, with distinct finite roots.',flush=True)
    if args.validate_only:
        print(f'Wrote {args.output}; production was not requested.',flush=True)
        return 0
    result['status']='production_running'
    for N in args.families:
        started=time.monotonic()
        DEADLINE=started+60*args.minutes_per_family
        ceiling=args.max_multiple or (10**9 if args.extend else (30 if N<=5 else 15))
        family=dict(N=N,k=2,data=[],status='running')
        result['families'].append(family)
        previous=None
        old={row['L']:row for group in saved.get('families',[]) if group['N']==N
             for row in group.get('data',[])}
        print(f'\nProduction N={N}, k=2',flush=True)
        for multiple in range(2,ceiling+1):
            L=N*multiple
            try:
                if L in old:
                    row=old[L]
                    z=np.array([complex(*v) for level in row['roots_by_level'] for v in level])
                    d=diagnostics(z,N,2,L)
                    rr=max(map(max,product_check(z,N,2,L)))
                    if (rr>=2e-9 or not np.isfinite(rr) or abs(d['E']-row['E'])>1e-10
                        or abs(d['energy_imag'])>=1e-9 or d['conjugation_error']>=1e-8
                        or d['reflection_error']>=1e-8 or d['minimum_separation']<=1e-7):
                        raise RuntimeError(f'Resume verification failed at N={N}, L={L}')
                else:
                    z,row=solve_ground(N,2,L,previous)
                row['A']=-6*L*(row['E']-bulk_energy(N,2)*L)/np.pi
                family['data'].append(row)
                previous=(L,z)
                print(f'L={L:4} E={row["E"]:.12f} A={row["A"]:.9f}'
                      f' residual={row["max_log_product_residual"]:.2e}'
                      f' time={row["elapsed_seconds"]:.2f}s',flush=True)
                if L not in old:
                    write_results(args.output,result)
            except (RuntimeError,TimeoutError) as exc:
                family.update(status='budget_exhausted' if isinstance(exc,TimeoutError) else 'solver_failed',
                              stopped_at_L=L,reason=str(exc))
                print(str(exc),flush=True)
                break
        else:
            family['status']='size_ceiling_reached'
        DEADLINE=float('inf')
        family['elapsed_seconds']=time.monotonic()-started
        family['fit']=fit_family(family['data'],N)
        write_results(args.output,result)
    result['families'].sort(key=lambda family:family['N'])
    result['status']='completed' if all(f['status']!='solver_failed' for f in result['families']) else 'partial_solver_failure'
    write_results(args.output,result)
    print_summary(result)
    print(f'Wrote all roots, energies, checks, and individual fits to {args.output}',flush=True)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
