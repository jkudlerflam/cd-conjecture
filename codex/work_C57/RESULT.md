# C57 result

[proved partial results; target unproved] Estimate (b) is not established here for interacting XXZ at fixed nonzero anisotropy, and therefore not for its retarded or multichannel extension. The work supplies a precise four-history bootstrap, a first-order XXZ thermal calculation, and a primary-source audit with a directly relevant prior-work lead. No unconditional thermal c=1 theorem is claimed.

## Literature finding

[read, announcement only] The 17 December 2024 Rutgers program lists Bruno Renzi, with Alessandro Giuliani, Vieri Mastropietro and Venanzio Alessio Santini, announcing a fermionic RG proof of the universal leading thermal free-energy correction for weakly interacting quantum chains. See [short talk B12, page 10](https://www.math.rutgers.edu/images/127_SMC_short_talks.pdf#page=10). The abstract gives no Hamiltonian, remainder estimate or theorem hypotheses. No full manuscript was located in the search. Its applicability to off-half-filled XXZ remains unverified. This is a useful lead and rules out saying that no such work has been announced.

[read, precise published analog] Giuliani–Mastropietro, [1211.3030](https://arxiv.org/pdf/1211.3030), equation (3.90) and the paragraph following (3.91), bound both history and geometric scalar differences by

    C |lambda| q_j^(2+theta/2) (1/(ell q_j))^(3-epsilon).

Together with (3.84) and (3.87), this is the shape of P with p>2. It is proved for their perturbed Ising model. Its free-Majorana reference and vanishing localized quartic cannot be replaced by the XXZ Luttinger reference without an additional argument.

[read, limits of other sources] BM 0907.2837, (81)–(82), controls the asymptotic beta function and running coupling; its Theorem 1.1 concerns exponents and susceptibility. BFM 1303.3681, Section 3.1, constructs the scalar expansion and its zero-temperature limit. MP 2108.13708, Section 8.2 and Theorem 8.5, provides analytic free energy and reference beta cancellation. None of these inspected estimates states P for XXZ. The recent rigorous XXZ quantum-transfer-matrix paper 2305.06679 explicitly retains conjectures in its physical free-energy identification, Theorem 1.4. Details and source links are in [SOURCE_AUDIT.md](SOURCE_AUDIT.md).

## What was proved

[proved conditional theorem] [PROOF.md, Section 3](PROOF.md) gives a bootstrap for the lattice/reference data X,Y at T and zero. The needed mixed difference is

    d_h=(X_h^T-Y_h^T)-(X_h^0-Y_h^0).

Uniform one-step derivative bounds and a zero-temperature defect O(q_h^theta) imply

    ||d_h|| <= C q_h^theta (T/q_h)^p.

The scalar comparison contains three pieces: the mixed state difference, a Hessian times the reference thermal shift times the zero-temperature defect, and the fixed-state thermal difference. The exact decomposition preserves both gains and yields P, including the final thermal shell, under explicitly stated hypotheses. A weighted-history version is also proved. The outstanding issue is establishing those uniform hypotheses for the renormalized XXZ expansion with all Gaussian scalar factors. A single-scale estimate with a field-weight loss cannot be iterated infinitely as stated.

[proved first perturbative coefficient] For hopping t and interaction V(n_x-1/2)(n_(x+1)-1/2), the thermodynamic first Taylor coefficient is

    f_1(T,mu)=(n_T-1/2)^2-|m_T|^2.

At fixed density rho=p/pi, direct Sommerfeld expansion gives

    f_1(T,rho)-f_1(0,rho)=T^2/(12t^2)+O(T^4).

This equals pi u_1 T^2/(6v^2), with v=2t sin p and first-order sound-speed shift u_1=2sin^2(p)/pi. Section 6 derives the distinct fixed-field coefficient and its agreement with the corresponding shifted sound speed. These are coefficientwise statements. There is no bound here on the higher interaction orders after division by T^2, so they do not prove the fixed-V thermal asymptotic.

## Retardation and multiple channels

[proved conditional extension] Exponentially small thermal matching errors can meet the initial thermal bounds. Static retarded endpoints must first rematch the field, forward coupling and frequency residue; their shifts can be linear in the endpoint size. C55's available progress still leaves the nonlinear all-scale estimates open. Its eventual action/source theorem would need the additional mixed scalar comparison used here.

[proved conditional extension] For finitely many channels the same bootstrap uses matrix/Banach norms. Its constants may depend on minimum speed, detuning and channel count. With a full comparison and reference normalization, C51 gives

    f(T)=e0-(pi/12)sum_a(1/u_(a,R)+1/u_(a,L))T^2+o(T^2).

The speeds are dressed Gaussian eigenmode speeds. Equality of right/left speeds needs a symmetry. The full physical-length ensemble transfer remains C58's separate task.

[numerical] `python3 -B checks.py` passes small quadrature and nonlinear-history controls. The initial finite-temperature threshold and floating-point cancellation failures were corrected and documented in PROOF Section 9. These checks test the partial mathematics, not an interacting RG flow. All files stay in work_C57; [STATE.md](STATE.md) is the local handoff.
