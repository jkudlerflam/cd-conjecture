"""Small finite algebra controls for C66. Run only from work_C66."""
from pathlib import Path
import numpy as np
from scipy.integrate import quad

out = ["# C66 finite checks", "", "[numerical] All checks below are finite algebra or normalization controls. They do not establish Lemma N, the comparison estimate, or a microscopic central charge.", ""]

def report(name, value):
    out.append(f"[numerical] {name}: {value}.")
    out.append("")

I = np.eye(2)
a = np.array([[0., 1.], [0., 0.]])
z = np.diag([1., -1.])
c0, c1 = np.kron(a, I), np.kron(z, a)
n0, n1 = c0.T @ c0, c1.T @ c1
Id = np.eye(4)

def bond(co):
    aa, bb, cc, uu, tr, ti = co
    tau = tr + 1j * ti
    return (aa * Id + bb * (n0 - Id/2) + cc * (n1 - Id/2)
            + uu * (n0-Id/2) @ (n1-Id/2)
            + tau * c0.T @ c1 + tau.conjugate() * c1.T @ c0)

coeff = np.sin(np.arange(1, 25).reshape(4, 6)) / 13
q = [bond(row) for row in coeff]  # 00,10,01,11
A, B, C, D = q[0], q[1]-q[0], q[2]-q[0], q[3]-q[1]-q[2]+q[0]
proj = [(Id-n0)@(Id-n1), n0@(Id-n1), (Id-n0)@n1, n0@n1]
lhs = sum(np.kron(pr, op) for pr, op in zip(proj, q))
rhs = np.kron(Id, A)+np.kron(n0, B)+np.kron(n1, C)+np.kron(n0@n1, D)
err = np.max(np.abs(lhs-rhs))
assert err < 1e-14
report("Four length-projector bond identity, maximum residual", f"{err:.3g}")
assert abs(np.trace(D)/4 - (coeff[3, 0]-coeff[1, 0]-coeff[2, 0]+coeff[0, 0])) < 1e-14

# Compare the analytic dictionary with derivatives of the independent
# Gaussian energy polynomial, for six unordered pairs and several seas.
coA = coeff[0]
coB = coeff[1]-coeff[0]
coC = coeff[2]-coeff[0]
coD = coeff[3]-coeff[1]-coeff[2]+coeff[0]

def ep(co, density, kap):
    aa, bb, cc, uu, tr, ti = co
    return aa+(bb+cc)*(density-.5)+uu*((density-.5)**2-abs(kap)**2)+2*np.real((tr+1j*ti)*kap)

def hp(co, density, kap, momentum):
    _, bb, cc, uu, tr, ti = co
    return bb+cc+2*uu*(density-.5)-2*uu*np.real(kap.conjugate()*np.exp(1j*momentum))+2*np.real((tr+1j*ti)*np.exp(1j*momentum))

maxerr = 0.
num = 0
for kp, kl, phi in [(1.1, .9, .2), (np.pi/2, np.pi/3, 0.), (1.3, 1.7, -.3)]:
    den, nu = kp/np.pi, kl/np.pi
    kap = np.exp(-1j*phi)*np.sin(kp)/np.pi
    kapl = np.sin(kl)/np.pi
    momenta = np.array([-phi+kp, -phi-kp, kl, -kl])
    def energy(x):
        npay = den+x[0]+x[1]
        nlen = nu+x[2]+x[3]
        cp = kap+np.dot(x[:2], np.exp(1j*momenta[:2]))
        cl = kapl+np.dot(x[2:], np.exp(1j*momenta[2:]))
        return ep(coA, npay, cp)+nlen*ep(coB+coC, npay, cp)+(nlen*nlen-abs(cl)**2)*ep(coD, npay, cp)
    predicted = {(0,1):4*np.sin(kp)**2*(coA[3]+nu*(coB[3]+coC[3])+(nu**2-abs(kapl)**2)*coD[3]),
                 (2,3):4*np.sin(kl)**2*ep(coD, den, kap)}
    for i in (0,1):
        for j in (2,3):
            hl = 2*nu-2*np.real(np.conjugate(kapl)*np.exp(1j*momenta[j]))
            predicted[i,j] = hp(coB+coC, den, kap, momenta[i])+hl*hp(coD, den, kap, momenta[i])
    step = .02
    for (i,j), pred in predicted.items():
        val = 0.
        for si in (-1,1):
            for sj in (-1,1):
                x = np.zeros(4)
                x[i], x[j] = si*step, sj*step
                val += si*sj*energy(x)
        actual = val/(4*step*step)
        maxerr = max(maxerr, abs(actual-pred))
        num += 1
assert maxerr < 1e-11
report(f"Six-coupling dictionary, {num} polynomial derivative comparisons, maximum residual", f"{maxerr:.3g}")

integral = quad(lambda t:t*t*np.exp(-8*t), 0, np.inf, epsabs=1e-13)[0]
assert abs(integral-1/256) < 1e-14
eta, mu, nu = .01, 1e-5, 1/3
phi = np.arctan(eta)
tp = np.sqrt(1+eta*eta)/256
kp = np.arccos(mu/(2*tp))
roots = [-phi+kp, -phi-kp]
slopes = [2/256*np.sin(p)+2*eta/256*np.cos(p) for p in roots]
vp = np.sqrt(4*tp*tp-mu*mu)
vl = 2*np.sin(np.pi*nu)
assert max(abs(slopes[0]-vp), abs(slopes[1]+vp)) < 1e-15
report("n=1 Gram integral residual", f"{abs(integral-1/256):.3g}")
report("Illustrative scheduled reference (eta=0.01, mu_p=1e-5, nu=1/3), v_p, v_l, v_l/v_p", f"{vp:.12g}, {vl:.12g}, {vl/vp:.12g}")

jl = 1.
for V in (-2., -2.3, -6.3):
    hb = -jl*(c0.T@c1+c1.T@c0)+V*n0@n1-V/2*(n0+n1)
    spec = np.linalg.eigvalsh(hb)
    assert np.min(spec) > -1e-14
    assert np.max(np.abs(spec-np.sort([0.,0.,-V/2-jl,-V/2+jl]))) < 1e-14
report("Attraction threshold shifted bond spectra at V=-2,-2.3,-6.3 with j_l=1", "all nonnegative, matching 0,0,-V/2-j_l,-V/2+j_l")

# Weighted Gaussian stability and absence of singularity at equal speeds.
mins = []
for vltest in (vl, vp, 2*vp):
    vabs = np.array([vp,vp,vltest,vltest])
    ag0 = np.diag(2*np.pi*vabs)
    root = np.diag(np.sqrt(2*np.pi*vabs))
    perturb = np.sin(np.arange(16).reshape(4,4)+.3)
    perturb = (perturb+perturb.T)/2
    perturb *= .2/np.linalg.norm(perturb,2)
    ag = ag0+root@perturb@root
    eig, vec = np.linalg.eigh(ag)
    assert np.min(eig)>0
    sq = (vec*np.sqrt(eig))@vec.T
    speeds = np.linalg.eigvalsh(sq@np.diag([1,-1,1,-1])@sq/(2*np.pi))
    assert np.sum(speeds>0)==2 and np.sum(speeds<0)==2
    mins.append(np.min(np.abs(speeds)))
report("Two positive and two negative Gaussian speeds, including coincident bare speeds", "passed in three weighted examples")

kappa, theta = .9, .3
const = 2**(-(kappa-theta))/((1-2**(-kappa))*(1-2**(-(kappa-theta))))
numsum = sum(2**(-kappa*d)*sum(2**(theta*r) for r in range(1,d+1)) for d in range(1,160))
assert abs(numsum-const)<1e-12
eta_rg, forcing, b, ae = .01, .002, 1.1, 1.7
assert b*eta_rg*const < .5
Bs = {}
for h in range(0,-101,-1):
    Bs[h] = ae*forcing*2**(theta*h)+b*eta_rg*sum(2**(-kappa*(k-h))*sum(Bs[j] for j in range(h+1,k+1)) for k in range(h+1,1))
ratio = max(Bs[h]/(2*ae*forcing*2**(theta*h)) for h in Bs)
assert ratio <= 1
report("C50 convolution constant residual and largest bound ratio over 101 scales", f"{abs(numsum-const):.3g}, {ratio:.9g}")

# Relevant-history recurrence: analytic power forcing with exact tail.
crel = 1/(2*(1-2**(-(1+theta))))
recerr = 0.
for h in range(-100,1):
    ah = -crel*forcing*2**(theta*h)
    ahm1 = -crel*forcing*2**(theta*(h-1))
    recerr = max(recerr, abs(ahm1-2*ah-forcing*2**(theta*h)))
assert recerr < 1e-17
report("Relevant-history geometric-tail recurrence residual", f"{recerr:.3g}")

pairerr = 0.
p = np.pi/3
for x in range(3,61,3):
    propagator = lambda y: np.sin(p*y)/(np.pi*y)
    corr = propagator(x)**2-propagator(x-1)*propagator(x+1)
    leading = 3/(4*np.pi**2*x*x)
    desc = 3/(4*np.pi**2*x*x*(x*x-1))
    pairerr = max(pairerr, abs(corr-leading-desc))
assert pairerr < 1e-15
report("Free pair descendant retained after Gaussian subtraction, maximum residual", f"{pairerr:.3g}")
out += ["[proved, scope] No external code or data was imported, and no interacting trajectory or microscopic c was computed.", ""]
Path('CHECKS.md').write_text('\n'.join(out))
print('\n'.join(out))
