"""Small analytic and quadrature checks for C95. Run only in this folder."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import math
from itertools import combinations
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit

HERE = Path(__file__).resolve().parent


def thermo(temp, tilt):
    def pressure(mu):
        return quad(lambda k: temp*np.logaddexp(0., (mu+2*math.cos(k))/temp),
                    0., math.pi, epsabs=2e-12)[0]/math.pi
    mu = brentq(lambda z: pressure(z)-tilt, -2., 2., xtol=1e-13)
    def integrand(k, power, response):
        f = expit((mu+2*math.cos(k))/temp)
        return math.cos(k)**power*(f*(1-f)/temp if response else f)
    def integral(power, response):
        return quad(lambda k: integrand(k, power, response), 0., math.pi,
                    epsabs=2e-12)[0]/math.pi
    rho, g1 = integral(0, False), integral(1, False)
    compress, g1_mu = integral(0, True), integral(1, True)
    a_rho = 1+g1*g1/rho**2-2*g1/rho*g1_mu/compress
    c = rho*compress*a_rho*a_rho/2
    return rho, c


def contact(n, v):
    size = n+v
    return n/size-math.sin(math.pi*n/size)**2/(n*size*math.sin(math.pi/size)**2)


def energy(n, v):
    size = n+v
    return -2*math.sin(math.pi*n/size)/math.sin(math.pi/size)


def toy_residual(coupling, temp):
    def logcosh(x):
        return np.logaddexp(x, -x)-math.log(2)
    delta = logcosh(math.sqrt(1+coupling**2/4)/temp)-logcosh(1/temp)
    return -temp*(np.logaddexp(0, -coupling/(2*temp)+delta)
                  -np.logaddexp(0, -coupling/(2*temp)))


def scalar_contact_coefficient(n, vacancies):
    size = n+vacancies
    shift = 0 if n % 2 else .5
    momenta = [2*math.pi*(j+shift)/size for j in range(size)]
    energies = [-2*math.cos(k) for k in momenta]
    filled = sorted(range(size), key=lambda j: energies[j])[:n]
    empty = set(range(size))-set(filled)
    coefficient = 0.
    for h1, h2 in combinations(filled, 2):
        for p1 in empty:
            p2 = (h1+h2-p1) % size
            if p2 not in empty or p2 <= p1:
                continue
            amplitude = 2/size*(math.cos(momenta[p1]-momenta[h1])
                                -math.cos(momenta[p2]-momenta[h1]))
            denominator = energies[p1]+energies[p2]-energies[h1]-energies[h2]
            coefficient += amplitude**2/denominator
    return coefficient/n


def main():
    c0 = 2/math.pi*(.25+1/math.pi**2)**2
    lines = ["# C95 checks", "", "[Numerical] Small quadratures and closed-form checks, with t=1. These do not test the full d=13 thermal trace.", "",
             f"Half-filling density-mode coefficient: {c0:.14g}.", "",
             "| T | s | Material density | Density-mode coefficient |",
             "| :-- | :-- | :-- | :-- |"]
    for temp in [.1, .03, .01, .003]:
        for tilt in [.5, 2/math.pi, .75]:
            rho, c = thermo(temp, tilt)
            assert 1/3 < rho < 2/3 and 0 < c < 4
            lines.append(f"| {temp:g} | {tilt:.9g} | {rho:.10g} | {c:.10g} |")
    assert abs(thermo(.003, 2/math.pi)[1]-c0) < 2e-5
    lines += ["", "[Numerical] At v=N, the adjacent-sector contact jump is O(1/N). The exact formulas are evaluated directly.", "",
              "| N | Crossing tilt | N times contact jump |",
              "| :-- | :-- | :-- |"]
    for n in [15, 30, 60, 120, 240, 480, 960]:
        jump = n*(contact(n,n)-contact(n,n+1))
        tilt = energy(n,n)-energy(n,n+1)
        assert abs(jump) < 1
        lines.append(f"| {n} | {tilt:.12g} | {jump:.12g} |")
    assert abs(jump-(.25+1/math.pi**2)) < .001
    lines += ["", "[Numerical] Exact two-particle-hole form-factor sum for the uniform contact. The N=4,6 and (N,v)=(9,4) values agree with the signed sums of C94's independently computed gap-resolvent rows.", "",
              "| N | Vacancies | Signed contact-kernel sum |",
              "| :-- | :-- | :-- |"]
    inherited = {(4,4): .02453690463133807,
                 (6,6): .024972535151530324,
                 (9,4): .008238644104655227}
    for n, v in [(3,3), (4,4), (6,6), (9,4), (15,15), (30,30)]:
        coeff = scalar_contact_coefficient(n,v)
        assert 0 <= coeff <= 48
        if (n,v) in inherited:
            assert abs(coeff-inherited[n,v]) < 2e-12
        lines.append(f"| {n} | {v} | {coeff:.14g} |")
    lines += ["", "[Numerical] Gapped two-level bath diagnostic, J=0.01. R2 is the Taylor coefficient at J=0 multiplied by J squared. The remainder is not uniformly smaller than R2 as T/J decreases.", "",
              "| T/J | R / J squared | (R-R2) / abs(R2) |",
              "| :-- | :-- | :-- |"]
    j = .01
    for ratio in [10, 1, .1, .01]:
        temp = ratio*j
        residual = toy_residual(j,temp)
        second = -math.tanh(1/temp)*j*j/16
        assert -j*j/8-1e-12 <= residual <= 1e-12
        lines.append(f"| {ratio:g} | {residual/j**2:.12g} | {(residual-second)/abs(second):.12g} |")
    assert abs((residual-second)/abs(second)-1) < 1e-8
    lines += ["", "[Numerical] All checks passed. The signed contact-sum bound, density-mode bounds, and two-level identities are proved separately in PROOF.md. No spatial first-moment bound or all-order d=13 estimate is inferred from these checks.", ""]
    (HERE/"CHECKS.md").write_text("\n".join(lines))
    print("All checks passed. See CHECKS.md.")


if __name__ == "__main__":
    main()
