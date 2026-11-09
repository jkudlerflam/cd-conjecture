"""Exact zero-energy resonance of the unlike pair at phi = pi/2, K = pi, when V3 = 0 and V2 = 2 t2 = -2.
det of the alpha-block of the boundary rows vanishes identically in u = |w| (so A has a pole), for general V2, V3."""
import sympy as sp
u=sp.symbols('u',positive=True); V2,V3=sp.symbols('V2 V3')
c=-u/4; s=sp.sqrt(1-c**2)
cosk=lambda n: sp.chebyshevt(n,c)
sink=lambda n: s*sp.chebyshevu(n-1,c) if n>=1 else 0
modes=[lambda n: cosk(n), lambda n: sink(n), lambda n: n*cosk(n), lambda n: n*sink(n)]
M=sp.Matrix([[ -2*u*m(1)-2*m(0)+V2*m(2) for m in modes],[ -2*m(1)+V3*m(3) for m in modes]])
dA=sp.factor(sp.simplify(M[:, :2].det())); dB=sp.factor(sp.simplify(M[:, 2:].det()))
print("det(alpha block) =", dA)
print("det(beta block)  =", dB)
print("at V2=-2, V3=0:", sp.simplify(dA.subs({V2:-2,V3:0})))
print("at V2=-2, V3=-2/3:", sp.factor(sp.simplify(dA.subs({V2:-2,V3:sp.Rational(-2,3)}))))
A=sp.simplify(M[:, :2].solve(M[:, 2:]))
ev=[sp.simplify(e) for e in A.eigenvals()]
print("eigenvalues of A (real basis) general V2,V3:"); [print("  ",e) for e in ev]
