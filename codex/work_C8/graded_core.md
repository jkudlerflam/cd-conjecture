# A uniform frustration-free fermion code, and its channel budget

Status: exact finite-dimensional construction and port formulas, proved
below. Its large-cell Clifford approximation is controlled. No finite-
coupling CFT or counterexample is asserted.

## 1. Repeated bond and complete ground code

Take q=n+1 local states, one even vacuum |0> and n odd colors |a>,
a=1,...,n. Define the graded swap and the positive bond

    P_g|a,b>=(-1)^(p(a)p(b))|b,a>,   h0=(I-P_g)/2,
    p(0)=0, p(a)=1 for a>0.                               (1)

On an open chain of N>=n sites, the complete zero space of the repeated
h0 has dimension 2^n. Label its orthonormal basis by subsets S of colors.
For r=|S|, its vector |S;N> is the signed equal superposition of every
configuration with precisely these r colors and N-r zeros. The sign is
the parity of the spatial permutation of the occupied colors relative
to their increasing order. Every amplitude has magnitude

    [N!/(N-r)!]^(-1/2).                                   (2)

Proof of completeness. Each adjacent graded swap must have eigenvalue
+1 in a zero vector. These swaps move any configuration with distinct
colors into any other configuration with the same color set. They fix
all relative coefficients to (2). A repeated odd color can be brought
to adjacent equal letters, where the graded swap is -1. Its coefficient
must vanish. There are no further constraints on the distinct-color
sets, and (2) explicitly satisfies every swap. Thus there is one ground
vector per subset. QED.

This is a genuinely uniform, nearest-neighbor, positive parent. Its
rank count includes n forbidden identical-color pairs, n vacuum-color
antisymmetric pairs, and n(n-1)/2 color-color symmetric pairs. Hence

    rank(h0)=n(n+3)/2.                                   (3)

## 2. Both endpoint maps

Identify the logical subset basis with the Fock space of n fermions.
Let f_a remove color a with the sign (-1)^(number of smaller colors),
let R=sum_a f_a^dagger f_a, and P=(-1)^R. Define

    g_N(R)=sqrt(N-R)/N,
    a_a=g_N(R) f_a,
    A_(a,x)=a_a+a_a^dagger,
    A_(a,y)=-i(a_a-a_a^dagger).                            (4)

The physical port operators are

    X_a=|0><a|+|a><0|,
    Y_a=-i|0><a|+i|a><0|.                                (5)

Their compressions at the first site are A_(a,x), A_(a,y).
To check the coefficient, remove a from a state with r occupied colors.
The surviving configurations with a fixed at site 1 number
(N-1)!/(N-r)!. Dividing by both normalizations from (2) gives

    [(N-1)!/(N-r)!] /
       sqrt([N!/(N-r)!][N!/(N-r+1)!])
      = sqrt(N-r+1)/N.                                  (6)

The fermionic sign is the number of colors smaller than a. This is
just the action of g_N(R)f_a, since R is evaluated after removal.
Creation follows by Hermitian conjugation.

At the last site the removal sign counts the colors greater than a,
so it differs by (-1)^(r-1). The two last-site compressions are

    B_(a,x)=-i A_(a,y) P,
    B_(a,y)= i A_(a,x) P.                                (7)

Consequently the physical inter-cell bond

    v=sum_a [Y_(a,right) tensor X_(a,left)
             -X_(a,right) tensor Y_(a,left)]              (8)

projects to sum_mu i A_mu P tensor A_mu, with the intended common
fermion-parity dressing. These are exact endpoint formulas, but for
n>1 the A_mu are not exact scaled Clifford matrices.

## 3. Controlled approximation and why it cannot violate c-d

Let Gamma_(a,x)=f_a+f_a^dagger and
Gamma_(a,y)=-i(f_a-f_a^dagger). For N>=n,

    ||sqrt(N) A_mu-Gamma_mu||
       = 1-sqrt(1-(n-1)/N).                              (9)

Each A_mu acts in independent 2 by 2 blocks joining a subset not
containing a to that subset with a added. The normalized amplitude in
such a block is sqrt(1-r_other/N), where 0<=r_other<=n-1. This proves
the equality in (9). For n=1 the compressions are exact at every N.

Thus this core realizes 2n approximate Clifford channels, with error
O(n/N), and both ports have the desired parity relation. But its free
first-order target has only c_target=n=q-1, before paying any selector
cost. The code dimension 2^n itself rules out more than 2n even Clifford
generators. It cannot access the proposed window above q-1.

Allowing p>=1 even colors and n odd colors in the same graded-swap
construction makes the code dimension grow with N. It does not repair
the exact-channel budget. The complete ground space contains the
constant product of any one even color. The boundary-Schmidt theorem
then bounds every even exact Clifford family by 2q-2 channels, with
q=p+n. Its normalized Gram obstruction also prevents arbitrarily
accurate larger families at fixed channel count on that entire code.

Removing all even colors instead gives the exterior power of C^q.
It has no zero states beyond N=q and so does not make arbitrarily long
cores. Selecting a proper subspace of any of these codes by new terms
is a different construction; neither this exclusion nor the exact
endpoint formulas establish what those terms can achieve.

## 4. Checks and the remaining issue

`graded_core_checks.py` constructs the physical code, repeated bond,
and endpoint matrices at (n,N)=(1,3),(2,3),(2,5),(3,4),(3,5). It verifies
the code residual, both ports, the error formula (9), and the physical
bond rank. At sizes with Hilbert dimension at most 256 it also checks
the complete ground multiplicity by dense diagonalization.

This is a useful local replacement for an abstract code in the vacancy
construction, but it only saturates the unblocked channel budget. It
therefore stalls as a counterexample route. A useful high-channel core
must avoid the product vectors of this family and have a larger logical
space, while retaining a repeated local parent and controlled ports.
