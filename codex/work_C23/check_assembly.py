"""Small finite algebra and hierarchy checks for C23; no microscopic simulation."""
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np


def filters(n):
    k = n + 1
    columns = []
    hs = []
    for j in range(n):
        q = [0] * (j + 1)
        for l in range(j + 1):
            c = math.comb(j + 2, l) * math.comb(j + 1, j - l)
            for r in range(j - l + 1):
                q[l + r] += c * math.comb(j - l, r) * (-1)**(j-l-r)
        x = [0] * k
        for r, c in enumerate(q):
            x[r] += c
            x[r+1] -= c
        columns.append(x)
        hs.append(Fraction(j+1, (2*j+4)*(j+3)))
    for a in range(n):
        assert sum(columns[a]) == 0
        for b in range(n):
            gram = sum(Fraction(columns[a][e]*columns[b][f], e+f+2)
                       for e in range(k) for f in range(k))
            assert gram == (hs[a] if a == b else 0)
    r2 = max(2 * sum(Fraction(columns[a][e]**2, 1)/hs[a]
                    for a in range(n)) for e in range(k))
    rho2 = min(1/(16*k*k*r2), Fraction(1, 8*k*(k+2)*n))
    return columns, hs, r2, rho2


def vertices(n):
    pairs = [(a,b) for a in range(-2,3) for b in range(-2,3)
             if (a-b)%2 == 0 and a*a+b*b <= 8]
    vectors = []
    for m in (-1,0,1):
        for ps in itertools.product(pairs, repeat=n):
            norm = 4*m*m + sum(a*a+b*b for a,b in ps)
            if norm not in (4,8) or sum(b for a,b in ps)%2:
                continue
            v = (-m,) + tuple(b for a,b in ps) + tuple(-a for a,b in ps)
            if next(x for x in v if x) > 0:
                vectors.append(v)
    return np.array(vectors, dtype=np.int64)


def sample(n, vs):
    count = len(vs)
    den = 1000003
    delta = 1/(4*count)
    rng = np.random.default_rng(2300+n)
    for attempt in range(100000):
        coords = np.concatenate(([rng.uniform(0.125,0.875)],
                 rng.random(n), rng.uniform(1/(32*n),0.5-1/(32*n),n)))
        nums = np.rint(coords*den).astype(np.int64)
        residues = (vs @ nums) % den
        dist = np.minimum(residues, den-residues)
        if int(dist.min()) * 4*count > den:
            assert np.all(np.sum(np.abs(vs),axis=1) <= 8)
            return nums, den, int(dist.min()), attempt+1
    raise RuntimeError('No sample')


def resources(n, b, rho2, delta):
    # Conditional resource example, NOT a compiled cell.
    # Hypotheses: core Delta=1, wrapper gap >= g, visibility unchanged,
    # ideal two-placement marker overlap 1/4 survives the parent.
    aa = lambda l: 4*math.sin(math.pi/(2*l))**2
    sel = min(aa(b-1)-2*aa(b)+aa(b+1),
              aa(b)-2*aa(b+1)+aa(b+2))
    eta = min((aa(b-1)-2*aa(b)+aa(b+1))/2,
              (aa(b)-2*aa(b+1)+aa(b+2))/3)
    head = 4*math.sin(3*math.pi/(2*(b+1)))*math.sin(math.pi/(2*(b+1)))
    helix = 1/((b+1)*(2*b+1))
    g = min(sel, head, helix)
    overlap = sum(2/math.sqrt(b*(b+1))*math.sin(math.pi*i/b)*
                  math.sin(math.pi*i/(b+1)) for i in range(1,b))
    beta = overlap**2/4 * math.cos(1/(2*b))**(2*b)
    B = 3.0
    j = float(rho2)
    vmin = 2*math.sin(math.pi/(16*n))
    eps = delta**2*vmin/(10000*(n+1))
    lam = eps/1000*min(g/B,g*g*j/B**3,beta*g/B,beta*math.sqrt(g/j))
    W = j*lam*lam
    u = W/beta
    e3 = 100*B**3*lam**3/g**2
    emix = 16*B*lam*u/g
    eu = 8*u*u/g
    assert e3/W <= eps/10*(1+1e-12)
    assert emix/W < eps
    assert eu/W < eps
    m=b-1
    C=(2*m+4)*2**m+m-4
    gamma=1/(8*C*C)
    D=C*(n+2)  # Explicit illustrative energy-event norm assumption.
    ereq=min(eps*W/100,eps*j*g*g/(100*B*B))
    Gamma=D+max(D*D/ereq,400*B*lam*D/(eps*W))
    Jc=Gamma/gamma
    return dict(b=b,selector_margin=sel,selector_eta=eta,head_gap=head,
                helix_gap_bound=helix,assumed_gap=g,head_overlap=overlap,
                assumed_beta=beta,clock_vertices=C,initialized_gap_per_Jc=gamma,
                assumed_energy_event_norm=D,eps=eps,lambda_value=lam,
                logical_scale=W,marker_u=u,length_bandwidth=4*W,
                cubic_bound=e3,mixed_bound=emix,marker_second_bound=eu,
                clock_schur_budget=ereq,Jc=Jc,clock_gap_Gamma=Gamma,
                leakage_amplitude_bound=(B*lam+u)/g,
                relevant_growth_budget=eps/delta**1.25,
                marginal_growth_budget=eps*math.log(1/delta))

out={}
for n in (2,3):
    cs, hs, r2, rho2=filters(n)
    vs=vertices(n)
    M=1+8*n+32*math.comb(n,2)+96*math.comb(n,3)+(128*math.comb(n,4) if n>=4 else 0)
    assert len(vs)==M
    nums,den,minnum,tries=sample(n,vs)
    coords=nums/den
    theta=2*math.pi*coords[1:n+1]
    widths=2*math.pi*coords[n+1:]
    out[str(n)]=dict(raw_columns=cs,h=[str(h) for h in hs],rmax_squared=str(r2),
        rho_squared=str(rho2),rho_squared_float=float(rho2),M=M,delta=1/(4*M),
        angular_margin=2*math.pi/(4*M),coordinate_shift_budget=1/(64*M),
        rational_coordinate_numerators=nums.tolist(),denominator=den,
        min_normalized_detuning=minnum/den,search_attempts=tries,
        centers=theta.tolist(),widths=widths.tolist(),mu=(2*np.cos(widths)).tolist(),
        velocities=(2*np.sin(widths)).tolist(),
        resources=resources(n,n+6,rho2,1/(4*M)))
# Equalization counterexample: scalar rescaling invariant cross ratio.
D0,D1=1,2
K=np.array([[1/(2*D0),1/(D0+D1)],[1/(D0+D1),1/(2*D1)]])
assert abs(np.linalg.det(K)-1/72)<1e-15
out['single_pole_counterexample']=dict(matrix=K.tolist(),determinant=float(np.linalg.det(K)))
Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
