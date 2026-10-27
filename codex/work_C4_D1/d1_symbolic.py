"""Symbolic checks for D1; no data inputs."""
import sympy as s
q=s.symbols('q',nonzero=True)
qn=lambda j: (q**j-q**(-j))/(q-q**(-1))
k=s.diag(q,1,1/q); ki=k.inv()
e=s.Matrix([[0,1,0],[0,0,1],[0,0,0]]); f=e.T
et=s.kronecker_product(e,k)+s.kronecker_product(ki,e)
ft=s.kronecker_product(f,k)+s.kronecker_product(ki,f)
ms=[2,1,0,1,0,-1,0,-1,-2]
cas=(q+1/q)*ft*et+s.diag(*[qn(m)*qn(m+1) for m in ms])
vals=[qn(j)*qn(j+1) for j in range(3)]
p=[]
for j in range(3):
    a=s.eye(9)
    for l in range(3):
        if l!=j: a=a*(cas-vals[l]*s.eye(9))/(vals[j]-vals[l])
    p.append(a.applyfunc(s.factor))
h=s.I*(q**4+1)/(q**4-1)*p[1]+s.I*(q**6-1)/(q**6+1)*p[0]
for i,j in [(2,4),(4,2),(4,6),(1,1),(2,2)]:
    v=s.factor(h[i,j]); anti=s.factor(v+v.subs(q,1/q))
    print('IK',i,j,'h=',v,'h-hdag=',anti)
c,t=s.symbols('c t',real=True)
E=s.zeros(4); E[1,1]=c-s.I*t;E[2,2]=c+s.I*t;E[1,2]=E[2,1]=-1
A=s.kronecker_product(E,s.eye(2));B=s.kronecker_product(s.eye(2),E)
H=2*c*A-A*B-B*A
paulis=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
for a in range(4):
 for b in range(4):
  for d in range(4):
   coeff=s.expand(s.trace(H*s.kronecker_product(paulis[a],paulis[b],paulis[d]))/8)
   if coeff!=0: print('TL','IXYZ'[a]+'IXYZ'[b]+'IXYZ'[d],s.factor(coeff))
