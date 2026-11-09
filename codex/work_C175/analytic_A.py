"""Closed-form threshold data at phi = pi/2 (gauge f(r) = i^r g(r) makes the K = pi problem real).
g-recurrence: 2g(r+2) + 2u g(r+1) + 2u g(r-1) + 2g(r-2) + V(r) g(r) = E g(r), E = -4 - u^2/2, u = |w|.
Free modes at threshold: double roots z = e^{+-i kappa}, cos kappa = -u/4.
Boundary rows (r = 2, 3): -2u g(1) - 2 g(0) + V2 g(2) = 0,  -2 g(1) + V3 g(3) = 0."""
import sympy as sp
u,r=sp.symbols('u r',positive=True)
V2,V3=sp.Rational(-2),sp.Rational(-2,3)
c=-u/4; s=sp.sqrt(1-c**2)
# real basis via Chebyshev-like: cos(kappa r), sin(kappa r) with cos kappa = c
def cosk(n): return sp.chebyshevt(n,c)
def sink(n): return s*sp.chebyshevu(n-1,c) if n>=1 else (0 if n==0 else -s*sp.chebyshevu(-n-1,c))
modes=[lambda n: cosk(n), lambda n: sink(n), lambda n: n*cosk(n), lambda n: n*sink(n)]
rows=[[sp.simplify(-2*u*m(1)-2*m(0)+V2*m(2)) for m in modes],[sp.simplify(-2*m(1)+V3*m(3)) for m in modes]]
M=sp.Matrix(rows)
Ma=M[:, :2]; Mb=M[:, 2:]
Areal=sp.simplify(Ma.solve(Mb))   # alpha = -Areal beta in the real basis (cos, sin)
print("A in (cos, sin) basis:"); sp.pprint(sp.simplify(Areal))
# convert to (e^{+i kappa r}, e^{-i kappa r}) channel basis: cos = (e+ + e-)/2, sin = (e+ - e-)/(2i)
T=sp.Matrix([[sp.Rational(1,2),sp.Rational(1,2)],[1/(2*sp.I),-1/(2*sp.I)]])  # rows: coeffs of cos,sin -> columns e+,e-
# a vector of (e+,e-) coefficients x maps to (cos,sin) coeffs y with x = P y, P = [[1/2,1/(2i)],[1/2,-1/(2i)]]
P=sp.Matrix([[sp.Rational(1,2),1/(2*sp.I)],[sp.Rational(1,2),-1/(2*sp.I)]])
Ach=sp.simplify(P*Areal*P.inv())
print("A in channel basis (e^{+i kappa r}, e^{-i kappa r}):")
for i in range(2):
    for j in range(2):
        print(i,j, sp.simplify(sp.expand_complex(Ach[i,j])))
print("eigenvalues:", [sp.simplify(e) for e in Ach.eigenvals()])
for uv in (0.5,1,2):
    print(uv, [complex(sp.N(Ach[i,j].subs(u,uv))) for i in range(2) for j in range(2)])
