# Endpoint data lose infinite radiative symplectic rank on a characteristic cone

## Abstract

Characteristic gluing and symplectic sewing ask different questions. We show this on the actual linearized Einstein null constraints on an outgoing Minkowski annulus. In a fixed odd tensor harmonic the full endpoint sphere-record map has three residual transport entries of rank two on compact histories. Its vertical space contains infinitely many mutually symplectically orthogonal nondegenerate planes, even after adjoining any fixed finite set of linear measurements. The same model admits a positive reconstruction: explicit local orthonormal coefficients determine the shear coordinate in the H¹₀ completion, recover the solved transports continuously, and converge in the Einstein symplectic pairing. At a finite subdivision, the trace and local coefficients assemble the original history and make the reduced flux additive. The result is linear and sector-specific. Vacuum Raychaudhuri yields a strict second-order area deficit under fixed initial area and expansion, showing why the vertical directions do not automatically integrate to nonlinear equal-cap families.

## 1. The actual characteristic endpoint map

Use four-dimensional vacuum general relativity with Newton constant G and the outgoing Minkowski double-null cone annulus of Aretakis, Czimek and Rodnianski (ACR), with affine radius 1≤r≤2, fixed unit lapse, fixed normal calibration and marked angular labels transported by the generator. Fix both cap embeddings and the initial sphere tangent data to zero. No final matching gauge operation is invoked. Let Y be a real degree ℓ≥2 spherical harmonic, X=*dY, and T a unit-L², symmetric odd tensor harmonic proportional to ∇₍A X_B₎. Its trace and double divergence vanish, and V=div T is nonzero with div V=0. This follows from the unit sphere commutator: div ∇₍A X_B₎=[2−ℓ(ℓ+1)]X/2. Write the fixed-determinant conformal metric perturbation as b(r)T, initially with b smooth and supported strictly in the open annulus.

The linearized equations are ACR Lemma 2.16, in the characteristic setup of §2.7 and §4. The radius perturbation has vanishing second derivative and zero initial value and derivative. Both it and the outgoing expansion perturbation therefore vanish. The intrinsic Gauss curvature perturbation is zero by double divergence. The outgoing shear and curvature are (r²/2)b′T and −(r²b′)′T/2 in these covariant angular components. With torsion, incoming tracefree shear and incoming curvature respectively aV,dT,zT, the ACR angular adjoint D₂*V=λT gives λ=[ℓ(ℓ+1)−2]/2>0: it is minus the trace-free symmetrized gradient. The remaining nontrivial radial equations and their quadratures are
```
 (r²a)′ = ½r²b′;          a(r) = (2r²)⁻¹ ∫₁ʳ s²b′(s) ds,
 r(d/r)′ = 2λa + ½rb′;    d(r)/r = b(r)/2 + 2λ∫₁ʳ a(s)/s ds,
 r(z/r)′ = 2λ(d/r²−a/r); z(r)/r = 2λ∫₁ʳ (d(s)/s³−a(s)/s²) ds.
```
The incoming scalar expansion is zero by its transport equation: the curvature, outgoing and lapse scalars and div V sources vanish. The incoming scalar lapse coefficient and its required transversal derivative vanish by the next scalar transports, again using div V=0. Outgoing lapse data are fixed. Thus after removing the identically zero rows of the full C² sphere data, and the automatically zero cap-local shear and outgoing curvature jets of compact profiles, the endpoint map is E(b)=(a(2),d(2),z(2)). If E=0, uniqueness of the homogeneous three-row ODE approached backward from the terminal cap makes the transported entries zero on its empty collar as well. Initial and final full sphere data then both vanish, including the low-charge projections. This does not require that all gauge-dependent gluing charges be finite dimensional on general angular data.

The three entries have **exact rank two**. Write M(r)=∫₁ʳ s b(s)ds. Integration by parts in the first equation and substitution into the next two give
```
 a=b/2−M/r²,   d/r=b/2+λM/r²,   z/r=2λ(λ+1)∫₁ʳ M(s)s⁻⁴ ds.
```
At the terminal cut b(2)=0, so a(2)=−M(2)/4 and d(2)/2=−λa(2). By interchanging the integration order, ∫₁² M(s)s⁻⁴ds = (1/3)∫₁² (t⁻²−t/8)b(t)dt. The two moments ∫₁² rb(r)dr and ∫₁² r⁻²b(r)dr are therefore exactly equivalent to imposing all three terminal transport entries. Their weight functions are nonproportional, so the two rows are independent on compact smooth profiles.

## 2. Einstein flux and its infinite-rank vertical defect

Use the Einstein–Hilbert covariant null potential, with the normalization and orientation for which the reduced spin-two flux is
```
 Ω(b,c) = (32πG)⁻¹ ∫₁² r² (bc′−cb′) dr.
```
Indeed the null volume density is r²dr dΩ; the conformally raised shear perturbation is b′T/2; the angular normalization is ∫|T|²dΩ=1. In the canonical decomposition of Hopfmüller and Freidel, the spin-one term is absent since the generator variation is fixed to zero and the spin-zero term is absent since the null volume and normal-scale variations are zero. The orientation of the cap normal alters the overall sign, not rank. The Einstein–Hilbert corner ambiguities relevant to this calculation vanish on the fixed gauge data and on the unperturbed cap collars of endpoint-vertical witnesses. This is a flux on reduced characteristic data; it does not identify all possible boundary potentials on non-fixed cuts. A nonzero compact b gives nonzero −(r²b′)′T/2 somewhere and hence cannot be a linearized diffeomorphism of Minkowski space.

**Theorem 1.** For every finite family of real linear measurement functionals M on the smooth compact radial histories, the restriction of Ω to ker(E,M) has infinite rank.

*Proof.* Set q=rb, so r²(bc′−cb′)=q p′−p q′ for p=rc. Choose a countable family of disjoint open radial subintervals with closures inside the annulus. In each interval select N disjoint nonzero compact smooth bumps F_j, with N>2(2+dim M). Stacking (E,M)(F_j/r) and (E,M)(F_j′/r) gives a real matrix after replacing E by its two independent moments, with fewer rows than columns. Choose a nonzero null vector; disjoint supports ensure F=Σc_jF_j is not zero. The two constrained histories F/r and F′/r both lie in ker(E,M). Their pairing is −(16πG)⁻¹∥F′∥²_L² by integration by parts. Histories attached to different intervals have zero cross pairings since the reduced flux is local in the shear coordinate even though solved transports are radial integrals. This supplies nondegenerate block-diagonal submatrices of order 2n for every n. The proof does not require continuity, independent rows, or locality of M. ∎

Any two-form pulled back from finite endpoint-plus-memory data has zero pairing on two vertical directions. The theorem therefore rules out that factorization, not existence of an endpoint-dependent section or a characteristic gluing solution. In particular choosing one control for every admissible pair of sphere records is not a presymplectic equivalence. The argument applies to the derivative of any smooth finite scalar-memory map defined at the background, but does not rule out a discontinuous encoding or an abstract infinite nonlinear compression.

## 3. Explicit local reconstruction and a commuting cut diagram

Define the reconstruction topology by q=rb∈H¹₀(1,2) with norm ∥q′∥₂. This is the completion of the smooth compact histories, not a claim that individual H¹ histories have classical C² sphere records. On an interval J=(u,v) of length h, write e_J for the unit-normalized Haar sign function, +h⁻¹/² on its left half and −h⁻¹/² on its right half, zero outside. Use all dyadic descendants of (1,2), including the interval itself. Put φ_J(r)=∫₁ʳ e_J(s)ds. These are compact tents in the closed interval J, zero on its boundary. Since the derivative is an isometry onto the mean-zero L² space, orthogonality and completeness of the Haar system on that space give the explicit local reconstruction
```
 c_J = ∫ q′ e_J,       q = Σ_J c_J φ_J,
 ∥q′∥₂² = Σ_J |c_J|².
```
Completeness also follows directly because dyadic step functions are dense in L², and their orthogonal decomposition uses the constant and finitely many Haar signs. Truncation at dyadic level m gives q_m with ∥q−q_m∥_{H¹₀}² equal to the omitted coefficient norm. The coefficients have a shear interpretation: the linear outgoing shear determines b′ on the marked generator, and b(1)=0 determines b and q; no metric reconstruction in the transverse spacetime is presumed.

Multiplication by 1/r is a bounded H¹ isomorphism on this annulus. The three transport quadratures define bounded linear maps H¹→H¹ for each solved coefficient, because their smooth weights and integration domains are bounded. In particular E extends continuously on the completion as the **transport** endpoint map. Smooth histories with zero transport rows are dense in its completed kernel: approximate q by compact smooth histories, then subtract fixed compact smooth representatives whose endpoint images form a basis of the finite-dimensional image. Such representatives exist because the dense image spans the same finite-dimensional subspace; their corrections vanish in norm with the residual rows. This is an option if exact row cancellation is desired at each sampling stage. H¹ convergence alone does not control b′ or b″ traces, so the phrase “full sphere record” remains reserved for the smooth compact histories. The tents are reconstruction functions, not new smooth classical vacuum solutions. Both Ω and E converge. In fact, writing Ω in q,p coordinates and using Cauchy–Schwarz and Poincaré yields
```
 |Ω(q,p)| ≤ (16πG)⁻¹ C_P ∥q′∥₂∥p′∥₂,
```
and the pairing error follows by inserting one tail at a time. There is no lower bound uniform over all frequency scales in this Hilbert norm.

Now cut at arbitrary marked radii 1=s₀<…<s_N=2. On each subinterval take the two H¹ traces of q and subtract their affine interpolant; the residual belongs to H¹₀ of that interval and admits the same explicit rescaled Haar reconstruction. Matched traces are precisely the condition for the assembled piecewise H¹ function to be global H¹: integration by parts in the weak derivative otherwise gives a delta proportional to each jump. Conversely the restrictions of a global H¹ function have equal traces. The squared derivative norm is the sum of the residual coefficient squares plus Σ_i |q(s_i)−q(s_{i−1})|²/(s_i−s_{i−1}), since the linear derivative is orthogonal to each residual's mean-zero derivative. Thus the reconstruction itself is continuous using that declared seam norm. Changing subdivision reconstructs the same q and p, and the outgoing shear flux integral splits exactly over subintervals. The solved triples are passed forward as their actual traces through the same triangular inhomogeneous ODE; zero-resetting them at intermediate cuts changes the solution.

Consequently the graph with restriction and assembly edges on the full trace-plus-coefficient data commutes and transports the reduced flux. Projecting it to only sphere records destroys the pairing by Theorem 1. This is a constructive sufficiency result for the stated topology and one harmonic, not a general statement that a nonlinear cut phase space is a Hilbert space, nor that one shear history reconstructs the spacetime off the null sheet. It does not quantify optimal bandwidth under higher regularity or inverse detector noise.

## 4. The first nonlinear consistency check

Take a smooth exact vacuum extension along the same marked affine generators, with the same **pointwise initial area density and expansion** through second order. At each fixed marked generator write $dA_\epsilon=\varphi_\epsilon^2d\Omega$, with $\varphi_\epsilon=r+\epsilon^2 w+O(\epsilon^3)$ because the first variation vanishes. The local affine vacuum equations $\theta=2\varphi_\epsilon^\prime/\varphi_\epsilon$ and $\theta^\prime=-\theta^2/2-|\sigma_\epsilon|^2$ give $\varphi_\epsilon^{\prime\prime}=-\varphi_\epsilon|\sigma_\epsilon|^2/2$, with $w(1,\omega)=w^\prime(1,\omega)=0$. The pointwise area factor is not the areal radius computed from the total area. The linearized mixed shear is b′T/2. Hence w″=−r(b′)²|T|²/8, and the coefficient of ε² in the physical area of the terminal sphere is
```
 −½ ∫₁² (2−r)r(b′)² dr
```
for the normalized harmonic. It is strictly negative for every nonzero compact shear. This is a statement conditional on an exact extension with the stated initial conditions, not an exact nonlinear endpoint reconstruction. The finite transport rows, the infinite linear flux, and the nonlinear area constraint therefore encode different information. On nonlinear characteristic data one must retain or solve the area constraint rather than matching two flat-cap histories solely because all first variations of their endpoint records agree.

## References and provenance

1. S. Aretakis, A. Czimek and I. Rodnianski, [*The characteristic gluing problem for the Einstein vacuum equations. Linear and non-linear analysis*](https://arxiv.org/abs/2107.02449), §2.7, Lemma 2.16 and §4. This supplies the characteristic constraints, not the symplectic rank obstruction.
2. S. Aretakis, A. Czimek and I. Rodnianski, [*The characteristic gluing problem for the Einstein equations and applications*](https://arxiv.org/abs/2107.02441). Overview of the separate gluing-existence problem.
3. A. Czimek and I. Rodnianski, [*Obstruction-free gluing for the Einstein equations*](https://arxiv.org/abs/2210.09663). Nonlinear gluing uses genuine controls and does not assert a flux descends through finite endpoints.
4. F. Hopfmüller and L. Freidel, [*Null conservation laws for gravity*](https://arxiv.org/abs/1802.06135), Eqs. (87)–(90). Null canonical potential used for the flux reduction.
5. M. P. Reisenberger, [*The symplectic 2-form for gravity in terms of free null initial data*](https://arxiv.org/abs/1211.3880). Closest free-data representation precedent; the local coefficient completeness here does not assert a new general null phase space.

The finite endpoint-vertical witness, the infinite-rank proof, and the local coefficient assembly refer to the same reduced sector and calibration. They do not identify a nonlinear cut algebra.
