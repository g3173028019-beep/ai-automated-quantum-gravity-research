# A boost obstruction for shell-flat spin-two memory infravacua

[Archive](../../README.md) · [Paper 23 in the index](../../docs/PAPER-INDEX.md#23-gravitational-memory-and-boosts) · [Research map](../../docs/RESEARCH-MAP.md)

**Archive connections:** [09](../../papers/09-canonical-einstein-scalar-hadamard/README.md) · [18](../../papers/18-positive-boundary-sewing/README.md). 09 concerns the existence of positive physical states and 18 separates positivity from extra compatibility requirements. 23 adds a concrete infinite-dimensional test: finite energy and positive quasifree data do not by themselves make a given symmetry implementable. Its null-infinity representation is not the state space of either earlier paper.

## Abstract

There is an explicit finite-Bondi-energy, ultraviolet-vacuum quasifree representation implementing the smooth electric-parity memory displacements in the following class, but no representation in the class implements even a single nontrivial standard boost if it implements a nonzero physical memory. The obstruction is a lower bound on the antilinear coefficient into *unsqueezed radial modes*, independent of angular ranks. The class fixes the radial subspace on each dyadic shell to the normalized flat profile in log frequency and allows arbitrary finite rotation-equivariant angular projectors and shellwise angular squeezing in the memory quadrature. We do not exclude squeezed smooth or overlapping radial bases, other quasifree representations, or dressed actions differing from the standard Lorentz automorphism.

## Physical radiative space and memory

On a chosen Bondi inertial frame at future null infinity take the two radiative helicities, with canonical one-particle norm proportional to \(\sum_{\lambda=\pm2}\int \omega\,d\omega\,d\Omega\,|f_\lambda(\omega,\Omega)|^2\); absorb the constant in \(f\). In terms of \(x=-\log(\omega/\omega_0)\) and \(w_\lambda=\omega f_\lambda\) this is the \(L^2(dx\,d\Omega)\) norm. Use real linear polarization tensors (the real STF angular tensor bundle) to fix the real structure. A memory step \(\Delta C_{AB}\) is the zero-frequency integral of Bondi news, and its corresponding displacement in this normalization has
\[
w(x,\Omega)=F(\Omega)+O(e^{-x}),\qquad x\to+\infty.
\]
An overall common phase is absorbed in defining the real memory quadrature. This follows either from the shear Fourier transform \(C(\omega)\sim \Delta C/(i\omega)\) with the radiative shear mode expansion \(C(\omega)\propto f(\omega)\), or from the leading \(1/\omega\) soft emission profile. Smooth finite-duration news, or an explicitly smooth ultraviolet completion of the displacement, ensures that the remainder is an ordinary one-particle vector.

For example in Bondi conventions \(\partial_u m=(D_A D_B N^{AB})/4-N_{AB}N^{AB}/8-4\pi G T_{uu}^{(2)}\). For an electric shear potential \(\Delta C_{AB}=(D_AD_B-\gamma_{AB}D^2/2)\Phi\), the left double derivative is \(D^2(D^2+2)\Phi/2\). The hard flux includes the gravitational news-squared term. Projecting out \(\ell=0,1\) and choosing stationary endpoint mass aspects with no higher-harmonic ordinary contribution gives a nonzero quadrupolar shear for hard flux differences with a quadrupolar component. This is physically accessible already with weak outgoing radiation: a smooth time-compact axisymmetric spin-two \(\ell=2,m=0\) linear news packet with zero time integral has squared norm proportional to \(\sin^4\theta\) and nonzero \(\ell=2\) projection. Its nonlinear order-squared electric memory is fixed by the integrated constraint. The \(\ell=0\) component is mass loss, not shear; reflection symmetry removes momentum recoil. Alternatively compare such a radiative balance with an isotropic hard matter energy flux of the same total energy. The result below needs just one such hard memory, does not require all formal smooth angular tensors to be realizable, and does not identify memory with an unconstrained extra bulk oscillator.

## Shell-flat squeeze and its energy

Write \(I_j=[jL,(j+1)L)\), \(L=\log 2\), for \(j\ge0\), and \(e_j=L^{-1/2}{\bf1}_{I_j}\). Let \(P\) be the orthogonal projection onto piecewise constant radial functions on these shells with unrestricted angular dependence, and zero at higher frequencies. The squeezed subspace on a shell is \(e_j\otimes V_j\), with \(V_j\) finite-dimensional, real and rotation invariant. On squeezed modes take the symplectic map
\[
 T(q+ip)=e^{-r_{j\ell}}q+i e^{r_{j\ell}}p
\]
in the memory-adapted real structure, allowing different strengths on electric/magnetic multiplets and, if desired, real orthogonal rotation-equivariant mixing within multiplicities. We use a diagonal common quadrature, \(r_{j\ell}\ge0\). \(T\) is the identity off chosen modes. The quasifree representation is Weyl\((Tf)\) on a Fock space over the transformed test completion. The dense test space may first use compact frequency support bounded away from zero. The displacement functional is implementable exactly when its transformed displacement \(Tw\) is in the completed one-particle space: continuity of the symplectic linear functional is necessary (test Weyl operators converge strongly at zero), and sufficient by its Riesz vector and Weyl translation. This distinction does not require a Fock vector for the unsqueezed \(w\).

For example take \(V_j\) to contain all tensor harmonics \(\ell\le j+2\) in both real tensor polarizations and \(r_j=\log(j+2)\). A smooth memory has transformed norm bounded in the soft part by a constant times
\[
\sum_j \left((j+2)^{-2}\|F\|^2+\|(1-P_{V_j})F\|^2\right)<\infty .
\]
The decaying \(O(e^{-x})\) remainder remains in the domain even if expanded in the enlarged quadrature, since \(e^{r_j}\) grows polynomially. The expected energy, defined as the sum of positive free radiative mode number times frequency relative to the ultraviolet vacuum, is bounded by
\[
\sum_j \dim V_j\,\omega_0 e^{-jL}\sinh^2 r_j<\infty .
\]
Indeed each one-mode covariance adds \(\sinh^2 r_j\) to its normal occupation, regardless of the anomalous two-point part. This gives finite energy as a covariance trace, not invariance of the squeezed vector under time translations. On compact radiative coordinate patches the covariance modification is smooth for the news and its derivatives: in log variables the news field insertion has an \(O(e^{-x})\) radial factor, time derivatives add powers of frequency, and a fixed number of angular derivatives of the finite harmonics costs only a polynomial in \(j\). The normal and anomalous covariance sums therefore converge after all fixed derivatives. No shear zero mode or individual graviton number finiteness is asserted.

## Doppler mixing and lower bound

For a fixed nonzero rapidity boost the action on on-shell waves is the massless helicity representation. If \(p'=\Lambda^{-1}p\), then \(\omega'=D(\Omega)\omega\) with \(D=\cosh\eta\pm\sinh\eta\cos\theta\) along the boost axis, and \(\Omega'=\alpha(\Omega)\) is aberration. In the \(w=\omega f\) normalization its action has the form
\[
 (U w)(x,\Omega)=D(\Omega)^{-1}R(\Omega)\,
 w(x-\log D(\Omega),\alpha(\Omega)).
\]
Here \(R\) is the real orthogonal rotation of the two linear graviton polarizations (diagonal phases of spin \(+2,-2\) in the helicity basis). The angular Jacobian is \(D^{-2}\), and the radial translation has unit Jacobian, so the action is unitary and real in the memory quadrature. In the equivalent forward-direction convention let K be the factor taking an incoming frequency to its boosted frequency, so K is the reciprocal of D evaluated at the boosted direction. The Fourier news and its zero-frequency memory then acquire one factor K and the spin rotation, while the integrated energy flux per solid angle acquires K³ (its total energy changes with the angular Jacobian K⁻²). The endpoint mass aspects must transform in the same Bondi balance: a rest-frame mass monopole is not a monopole in all frames. The quadrupole example is specified in a single inertial frame; it does not suppress these endpoint terms after a boost. The angle dependence, not a uniform energy rescaling, is essential. At sufficiently soft input shells their translated support stays away from the hard endpoint \(x=0\).

For a test input \(e_jY\) with a real unit spin-two harmonic \(Y\) independent of \(j\), the outgoing radial shape at each angle is the translated indicator, split between two adjacent shells (after an angle-dependent integral number of whole shells). If \(\tau(\Omega)=\{\log D(\Omega)/L\}\), the projection of this unit radial vector onto the entire piecewise constant radial subspace has squared norm \((1-\tau)^2+\tau^2\). Consequently
\[
 \|(1-P)U(e_jY)\|^2 =
 c_{\Lambda,Y}:=\int d\Omega\,D^{-2}|Y(\alpha(\Omega))|^2\,2\tau(1-\tau)>0.
\]
Integer shifts occur only on finitely many angular latitudes of zero measure for a fixed finite boost. The actual spin rotations cancel in the norm. This leakage probes the unsqueezed space even if every angular mode could otherwise be retained; increasing angular ranks cannot reduce it.

The conjugated action on the completed Fock space would have to be \(TUT^{-1}\) and satisfy the bosonic Shale condition (bounded canonical map, Hilbert–Schmidt antilinear part). For an included squeezed eigenmode \(g_j=e_jY\), both real quadratures map by \(U\), while \(T^{-1}\) multiplies them by \(e^{+r_j}\) and \(e^{-r_j}\) respectively. On output \(1-P\), \(T\) is the identity. The antilinear coefficient on \(g_j\), projected onto this complement, is thus exactly \(\sinh r_j (1-P)Ug_j\). Summing over these orthonormal input modes gives a lower bound \(c_{\Lambda,Y}\sum_j\sinh^2 r_j\) on the squared Hilbert–Schmidt norm; cancellations between outputs of different inputs cannot affect this norm. If arbitrary rotation-equivariant eigenvectors within the two polarizations are allowed, a fixed \(\ell\) has finite dimension and the positive leakage form above has a strictly positive minimum on its unit sphere. The same lower bound applies uniformly to its varying normalized real eigenvectors.

A physical memory with nonzero \(\ell=2\) component requires \(\sum_j\|T_j F_{\ell=2}\|^2<\infty\) in its real quadrature, since its shell coefficients approach \(\sqrt{L} F_{\ell=2}\) exponentially. On every sufficiently soft shell there must therefore be an eigenvalue \(e^{-r_j}\to0\) in that fixed angular block; absent or bounded squeezing on an infinite subsequence forbids implementability. Select one corresponding unit real eigenvector per shell. The boost Hilbert–Schmidt lower bound diverges. In fact the norms of the projected images of the unit real inputs under the full conjugated boost grow without bound, independently excluding a bounded implementation. Compact shell inputs suffice to establish that obstruction, without assuming the conjugated boost has any extension. A representation cannot possess the requested strongly continuous full proper Lorentz action if this individual boost is not implementable.

## Scope and comparison

The result identifies what fails: the Doppler shift takes the heavily contracted memory quadrature into radial shapes deliberately left at their ordinary vacuum variance. It is not ordinary Fock orthogonality, since the same class explicitly implements the physical memory at finite mean energy when boosts are not required. It is a no-go in the shell-flat radial construction, not for more general finite radial subspaces, nonquasifree representations, or a different (dressed) Lorentz action. Only the leading soft phase space and the physical memory displacement determined by the nonlinear flux balance enter; no all-orders scattering operator is constructed.

[Prabhu–Satishchandran–Wald](https://arxiv.org/abs/2203.14334) discuss gravitational memory representations and infrared scattering; [Kraus–Polley–Reents](https://www.numdam.org/item/AIHPA_1977__26_2_109_0.pdf) supply the electromagnetic infrared representation precedent. The specific computed claim here is the fixed spin-two shell leakage under the geometric boost, combined with the physical memory implementability condition and the finite-energy counterexample when that symmetry requirement is dropped.

[Prabhu and Satishchandran, *Infrared finite scattering theory: scattering states and representations of the BMS group*](https://doi.org/10.1007/JHEP08(2024)055) construct direct-integral representations over Lorentz orbits of memory sectors, with generators in that larger representation even when a fixed sector lacks them. This is not an implementation within one shell-flat quasifree representation, and our obstruction does not exclude their construction. The present distinction is an explicit Doppler leakage bound into the untouched radial complement despite finite energy and implementable shifts in a single shell-flat representation.
