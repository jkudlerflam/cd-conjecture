"""Symbolic checks for C150: kernel identities and the derivative recursions."""
import sympy as sp
s, x, eta, y = sp.symbols('s x eta y', real=True)
T = sp.sin(2*eta)/(sp.pi*(sp.cosh(2*s)-sp.cos(2*eta)))
Torig = sp.tan(eta)*sp.sech(s)**2/(sp.pi*(sp.tanh(s)**2+sp.tan(eta)**2))
theta = -2*sp.atan(sp.tanh(s)/sp.tan(eta))
a1 = sp.sin(eta)/(sp.pi*(sp.cosh(2*x)+sp.cos(eta)))
p = 2*sp.atan(sp.tanh(x)*sp.tan(eta/2))
e = -sp.sin(eta)**2/(2*(sp.sinh(x)**2+sp.cos(eta/2)**2))
import random
for _ in range(5):
    vals = {s: random.uniform(-2, 2), eta: random.uniform(0.01, 0.7), x: random.uniform(-2, 2)}
    print('T forms', float((T-Torig).subs(vals)), ' -theta\'/2pi - T', float((-sp.diff(theta, s)/(2*sp.pi)-T).subs(vals)),
          ' p\'/2pi - a1', float((sp.diff(p, x)/(2*sp.pi)-a1).subs(vals)), ' e + pi sin(eta) a1', float((e+sp.pi*sp.sin(eta)*a1).subs(vals)))
# ratios T^(i)/T as functions of C = cosh 2s
C, c = sp.symbols('C c', real=True)
for i in (1, 2):
    r = sp.simplify(sp.diff(T, s, i)/T)
    print('T^(%d)/T =' % i, r)
# derivative recursion for P_m = D^(m-1) [T(x-y)/sigma(y)], D = (1/sigma) d/dy
sig = sp.Function('sigma')(y)
Tf = sp.Function('T')
P = Tf(x-y)/sig
out = [P]
for m in range(2, 5):
    P = sp.expand(sp.diff(P, y)/sig)
    out.append(P)
for m, P in enumerate(out, 1):
    print('P_%d =' % m, P)
# same for f = e(y)
ef = sp.Function('e')(y)
F = ef
for m in range(1, 4):
    F = sp.expand(sp.diff(F, y)/sig)
    print('D^%d e =' % m, F)
