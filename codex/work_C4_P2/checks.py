from fractions import Fraction as F

def solve(a):
    n=len(a)
    m=[[F(v) for v in row]+[F(1)] for row in a]
    for i in range(n):
        p=next(j for j in range(i,n) if m[j][i])
        m[i],m[p]=m[p],m[i]
        z=m[i][i]; m[i]=[x/z for x in m[i]]
        for j in range(n):
            if j!=i:
                z=m[j][i];m[j]=[x-z*y for x,y in zip(m[j],m[i])]
    return [2*row[-1] for row in m]

def path(n):
    a=[[2*int(i==j)-int(abs(i-j)==1) for j in range(n)] for i in range(n)]
    return a

def arms(lengths):
    n=1+sum(lengths); a=[[2*int(i==j) for j in range(n)] for i in range(n)]
    nxt=1
    for length in lengths:
        prev=0
        for _ in range(length):
            a[prev][nxt]=a[nxt][prev]=-1;prev=nxt;nxt+=1
    return a

if __name__=='__main__':
    for r in range(2,13):
        a=path(r);assert min(solve(a))==r
        b=path(r);b[-2][-1]=-2
        assert min(solve(b))==min(2*r,r*(r+1)//2)
        c=path(r);c[-1][-2]=-2
        assert min(solve(c))==2*r-1
        if r>=4:assert min(solve(arms([1,1,r-3])))==2*r-2
    for r,branch,expected in [(6,[2,2,1],16),(7,[2,3,1],27),(8,[2,4,1],58)]:
        values=solve(arms(branch));assert min(values)==expected
        print('E'+str(r), 'C_i =',list(map(str,values)))
    f=path(4);f[1][2]=-2
    for name,a,expected in [('F4',f,16),('G2',[[2,-1],[-3,2]],6)]:
        values=solve(a);assert min(values)==expected
        print(name,'C_i =',list(map(str,values)))
    for n in range(2,30):
        for k in range(1,30):
            c=F(k*(n*n-1),k+n)
            assert c<=k*(n-1)
        assert F(2*n*(2*n+1),n+3)>2*n-1 # Sp(2n)_2
    print('All exact rational checks passed.')

# Check the exact one-defect dispersion of a mixed graded permutation.
# The background has odd parity. The defect has even parity.
if __name__=='__main__':
    import cmath
    import math
    for length in range(4,10):
        for mode in range(length):
            q=2*math.pi*mode/length
            v=[cmath.exp(1j*q*j) for j in range(length)]
            # (sum P_g + L I) restricted to this subspace.
            hv=[2*v[j]+v[(j-1)%length]+v[(j+1)%length] for j in range(length)]
            energy=2+2*math.cos(q)
            assert max(abs(a-energy*b) for a,b in zip(hv,v))<1e-12
    print('Mixed graded-permutation one-defect dispersion checks passed.')
