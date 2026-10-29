"""Small normalization checks for C109. Run only in this folder."""
from pathlib import Path
from fractions import Fraction
import math
import numpy as np
from scipy.integrate import quad

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
P = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)
report = ['# C109 checks', '', '[Numerical] All calculations are small deterministic normalization checks. The largest diagonalization is 64 by 64. They do not establish infrared matching.', '']

def record(s):
    report.append('[Numerical] PASS. '+s)
    report.append('')

def spread(h):
    e = np.linalg.eigvalsh(h)
    return (e[-1]-e[0])/2

def hxy(gamma):
    return -(1+gamma)*np.kron(X,X)/2-(1-gamma)*np.kron(Y,Y)/2-(np.kron(Z,I)+np.kron(I,Z))/2

def embed(a, i, b, j, n):
    ops = [I for _ in range(n)]
    ops[i], ops[j] = a, b
    out = np.array([[1]], dtype=complex)
    for op in ops:
        out = np.kron(out,op)
    return out

def chain(terms, edges, n):
    out = np.zeros((2**n,2**n), dtype=complex)
    for i,j in edges:
        for c,a,b in terms:
            out += c*embed(a,i,b,j,n)
    return out

def termsxy(gamma):
    return [(-(1+gamma)/2,X,X), (-(1-gamma)/2,Y,Y), (-.5,Z,I),(-.5,I,Z)]

def dispersion(k,gamma):
    return 2*np.sqrt((2*np.sin(k/2)**2)**2+gamma**2*np.sin(k)**2)

worst_spectrum = worst_gauge = worst_slope = 0.
for gamma in [1., .7, .2, .01, .001]:
    h = hxy(gamma)
    expected = np.array([-math.sqrt(1+gamma**2), -1., 1., math.sqrt(1+gamma**2)])
    err = np.max(np.abs(np.linalg.eigvalsh(h)-expected))
    assert err < 1e-12
    worst_spectrum = max(worst_spectrum,err)
    assert np.linalg.norm(P@h@P-h)<1e-12
    for o in [X, .3*Y+.7*Z, -.4*X+.2*Y+.1*Z]:
        b = np.kron(o,I)-np.kron(I,o)
        assert np.linalg.norm((h+b+P@(h+b)@P)/2-h)<1e-12
        diff = spread(h+b)-spread(h)
        worst_gauge = min(worst_gauge,diff)
        assert diff > -1e-12
    k = gamma*1e-4
    err = abs(dispersion(k,gamma)/(2*gamma*k)-1)
    worst_slope = max(worst_slope,err)
    assert err < 3e-9
record(f'XY two-site spectrum and swap-gauge identity, maximum eigenvalue error {worst_spectrum:.3e}; no tested gauge lowered the spread. Small-k slope relative error {worst_slope:.3e}.')

worst_ring=0.
for n in [4,6]:
    edges=[(i,(i+1)%n) for i in range(n)]
    for gamma in [.2,.7,1.]:
        H=chain(termsxy(gamma),edges,n)
        exact=np.linalg.eigvalsh(H)[0]
        k=(2*np.arange(n)+1)*np.pi/n
        jw=-np.sum(dispersion(k,gamma))/2
        err=abs(exact-jw)
        worst_ring=max(worst_ring,err)
        assert err<1e-11
record(f'Periodic spin-ring ground energies versus the independently evaluated Jordan-Wigner vacuum sum at L=4,6 and three anisotropies, maximum error {worst_ring:.3e}.')

worst_twist=0.
for n in [4,6]:
    H=chain(termsxy(1.),[(i,(i+1)%n) for i in range(n)],n)
    # Z seam reverses the one crossing XX interaction.
    Ht=H+2*embed(X,n-1,X,0,n)
    delta=np.linalg.eigvalsh(Ht)[0]-np.linalg.eigvalsh(H)[0]
    predicted=2*math.tan(math.pi/(4*n))
    err=abs(delta-predicted)
    worst_twist=max(worst_twist,err)
    assert err<1e-11
    assert np.linalg.norm(Ht-H,2)<=2*math.sqrt(2)+1e-12
record(f'Ising Z2 seam energy delta E=2 tan(pi/(4L)), maximum error {worst_twist:.3e}. This gives lim L delta E=pi/2, consistent with v=2 and x=1/8.')

terms=[(.3,X,X),(.2,Y,Y),(.17,Z,Z),(.11,X,Z),(.07,Z,I),(-.13,I,Y)]
ncopy,L=2,3
n=ncopy*L
bulk=[(i*ncopy+s,(i+1)*ncopy+s) for i in range(L-1) for s in range(ncopy)]
untwisted=bulk+[((L-1)*ncopy+s,s) for s in range(ncopy)]
twisted=bulk+[((L-1)*ncopy+s,(s+1)%ncopy) for s in range(ncopy)]
Hu=chain(terms,untwisted,n)
Ht=chain(terms,twisted,n)
Hring=chain(terms,[(i,(i+1)%n) for i in range(n)],n)
Hsmall=chain(terms,[(i,(i+1)%L) for i in range(L)],L)
specerr=np.max(np.abs(np.linalg.eigvalsh(Ht)-np.linalg.eigvalsh(Hring)))
energyerr=abs(np.linalg.eigvalsh(Hu)[0]-ncopy*np.linalg.eigvalsh(Hsmall)[0])
assert specerr<1e-11 and energyerr<1e-11
hb=sum(c*np.kron(a,b) for c,a,b in terms)
# This uses the unoptimized norm, sufficient to test the seam inequality itself.
assert np.linalg.norm(Ht-Hu,2)<=2*ncopy*np.linalg.norm(hb,2)+1e-11
record(f'Cyclic replica seam versus a single long ring, entire spectrum maximum error {specerr:.3e}; untwisted product ground-energy error {energyerr:.3e}. A generic complex, oriented bond was used.')

for N in range(2,201):
    c=Fraction(2*(N-1),N+2)
    dims=[Fraction(k*(N-k),N*(N+2)) for k in range(N)]
    avg=sum(dims)/N
    assert c/avg==Fraction(12*N,N+1)
    assert c/dims[1]==2*N
    assert c/dims[N//2]<=8
assert Fraction(1,2)/Fraction(1,8)==4
assert Fraction(4,5)/Fraction(2,15)==6
record('Parafermion twist identities verified using rational arithmetic for N=2 through 200, including the maximal-twist factor 8 and averaged factor below 12.')

worst_fourier=0.
for R in [.5,2.]:
    xs=np.arange(8.)
    kernel=1/(1+((xs[:,None]-xs[None,:])/R)**2)**2
    assert np.linalg.eigvalsh(kernel)[0]>-1e-12
    for q in [0.,.3,1.]:
        if q==0.:
            val=2*quad(lambda x:1/(1+(x/R)**2)**2,0,np.inf,epsabs=1e-11)[0]
        else:
            val=2*quad(lambda x:1/(1+(x/R)**2)**2,0,np.inf,weight='cos',wvar=q,epsabs=1e-11)[0]
        expected=np.pi*R/2*(1+R*abs(q))*np.exp(-R*abs(q))
        err=abs(val-expected)
        worst_fourier=max(worst_fourier,err)
        assert err<1e-9
record(f'Positive covariance countermodel: finite Gram matrices positive and Fourier-transform integration agrees with the analytic formula to {worst_fourier:.3e}.')

report.extend(['[Proved, analytic conclusions checked above] The tests neither prove a uniform seam budget nor a general bound on A_neutral. The velocity-floor counterexample and replica identity are established analytically in PROOF.md.', ''])
Path('CHECKS.md').write_text('\n'.join(report))
print('\n'.join(report))
