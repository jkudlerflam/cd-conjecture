# Interval-arithmetic (mpmath iv, 30 digits) version of pos_rigorous.py, added by the coordinating Claude session. Output: worst cell 0.938945 < 1, direct margin 0.0210 > 0, f_small(1) = 0.8041.
from mpmath import iv, mpf
iv.dps=30
coth=lambda x:(iv.exp(2*x)+1)/(iv.exp(2*x)-1)
c=2*iv.pi*iv.sin(iv.pi/8)
def br2(l):
    t=iv.exp(-2*l); return iv.mpf(1)/2+2*(t/(1+t)+t**2/(1+t**2)+t**3/(1+t**3))
worst=None
for k in range(36):
    l0=iv.mpf(1)+iv.mpf('0.36')*k/36; l1=iv.mpf(1)+iv.mpf('0.36')*(k+1)/36
    v=l1*(2*coth(l0)-br2(l1))/(c*(1-2*iv.exp(-iv.pi**2/l1))**2)
    worst=v if worst is None or v.b>worst.b else worst
print("worst cell upper enclosure:",worst.b)
l=iv.mpf('1.36'); print("direct margin enclosure:",((1-2/(iv.exp(l)+1))-coth(l)/2).a)
l=iv.mpf(1); print("f_small(1) upper:",((2*l*coth(l)-iv.log(2))/(c*(1-2*iv.exp(-iv.pi**2/l))**2)).b)
