# Infinite mass-coupling susceptibility on the neutral de Sitter static-patch algebra


## Abstract

Comparing massive Euclidean states with the massless current state on a fixed de Sitter static-patch algebra removes the homogeneous position pole, while retaining the small homogeneous momentum fluctuations. The resulting states are normal, yet their Araki relative entropy has infinite second susceptibility in the parameter $m^2$ at zero. The lower bound comes from bounded approximations to neutral local squares: their linear mass response grows logarithmically as the aperture narrows, as does their Kubo–Mori covariance. This result does not assert infinite relative entropy at any fixed nonzero mass.

## Statement

In four-dimensional de Sitter take the Euclidean states of the minimally coupled free real scalar, with squared mass parameter μ=m². Identify them at the reflection-symmetric global Cauchy equator by fixed canonical neutral Cauchy observables, and use the von Neumann algebra M of the massless current in the northern static patch (the domain of dependence of the open northern hemisphere). Then the massive neutral restrictions are normal on M, and converge on every fixed finite set of Cauchy Weyl probes to the massless Euclidean current state. Nevertheless, for Araki relative entropy,
```
liminf_{μ↓0}  2 S_M(ω_μ || ω_0)/μ² = +∞.                               (1)
```
This is an ultraviolet modular-information obstruction for the mass **coupling** μ, not a residual homogeneous momentum singularity. The lower bound uses only compactly interior, neutral field probes, but relies on the static-patch modular group of the fixed algebra M. It is not a proof about all smaller regions, every fixed μ>0, a different comparison Cauchy surface, or the differently parametrized Hessian in m.

## Common algebra, pole and normality

Let R be the de Sitter radius, Σ=S³_R the t=0 equator, and U its northern hemisphere. Use the induced volume and let π be the normal derivative at t=0. Weyl labels are (f,g) in C_c^∞(U)² with ∫Σ f=0, and pairing φ(f)+π(g) and its canonical symplectic form are identical for all μ. The massless observables are the time-slice representatives of spacetime tests invariant under φ↦φ+c: the position test is neutral, while the derivative test g is arbitrary. Equivalently, in the massless solution-space presentation the symplectic pairing against the constant solution vanishes; common labels at one Cauchy surface do not identify the dynamics of different masses. There is no local null direction: a compactly supported g orthogonal to all neutral f on U would be constant in U and hence vanish. On the whole compact Σ the corresponding constant g is instead a central momentum in this invariant CCR; the massless Euclidean current representation selects zero for it.

A useful Euclidean verification uses the Dirichlet-to-Neumann map K_μ of −Δ_{S⁴_R}+μ on a Euclidean hemisphere bounded by Σ, with outward derivative convention. The reflected Gaussian equatorial state has φ covariance (2K_μ)⁻¹ and π covariance K_μ/2, with vanishing symmetric cross term. This follows harmonic by harmonic by the Euclidean on-shell Gaussian ground-state wavefunctional exp(−(φ,K_μφ)/2). The constant harmonic normalized on Σ has k_0(0)=0 and
```
k_0(μ)=μ Vol(hemisphere of S⁴_R)/Vol(Σ)+O(μ²);
```
integrate the constant-harmonic elliptic equation in the hemisphere and use that the normalized solution is constant at μ=0. Thus its momentum variance is O(μ), not inverse-power divergent. Its position variance produces exactly the rank-one pole 1/(μ Vol(S⁴_R)) in the Euclidean covariance kernel. The Dirichlet problem itself is analytic at μ=0 (Dirichlet Laplacian has a gap); only inversion of the constant boundary eigenvalue is singular.

On the orthogonal complement of the constant mode the massless and massive equatorial pure Gaussians are Fock-equivalent. Indeed K_μ has order 1 with μ-independent principal and subprincipal symbols; K_μ−K_0 has order −1 with Cauchy parameter small. This follows by solving the elliptic collar factorization (the added order-zero bulk potential μ first contributes μ/(2|k|) to the normal symbol). Since spherical symmetry diagonalizes K_μ, k_ℓ(μ)/k_ℓ(0)−1=O(μ(1+ℓ)⁻²) for ℓ≥1, and Bogoliubov coefficients with the same decay are square-summable with multiplicities (ℓ+1)². This proves the required normality of the nonzero modes by the implementability condition. For the constant-momentum component of the massive neutral state, independently Gaussian with O(μ) variance, choose a smooth mean-zero function on Σ which is constant on the closed northern hemisphere and compensated strictly in the southern interior. Its position Weyl operator implements on the local momentum smearings precisely the constant shift, with appropriate normalization. The massive local neutral state is the Gaussian mixture of these unitarily displaced normal nonzero-mode states (global unitaries need not belong to M). This supplies normality on the *fixed* massless local von Neumann algebra, without assigning zero to the massive homogeneous momentum or declaring it a separately measurable local sector.

At the Euclidean two-point level, if L=−Δ_{S⁴_R} and P_⊥ projects off constants,
```
G_μ = (μ Vol(S⁴_R))⁻¹ + G_0^⊥ − μ (L|_⊥)⁻² P_⊥ + O(μ²),             (2)
```
where the first term denotes its constant kernel. In particular every fixed smooth neutral equator covariance is one-sided analytic in μ near zero; the same is true of any finite normal Gaussian probability measurement of these filtered fields. On the full field algebra even a single nonneutral position smear has variance growing as μ⁻¹, so there is no regular Gaussian massless limit there: the limiting Weyl characteristic function is discontinuous at zero. Its Gaussian measurement relative entropy against any fixed finite, positive reference variance diverges already at that one observable. This elementary full-field IR obstruction is distinct from (1).

## Modular test and its quantitative UV estimates

The massless current restriction ω_0 is faithful KMS for the static Killing time, inverse temperature β=2πR relative to central proper time. In the massive Euclidean theory this follows from the de Sitter wedge KMS property; for the massless current it follows on the gauge-invariant field distributions after the constant Euclidean term cancels, by the same regular angular continuation on S⁴. In a static tube strictly inside the horizon the Euclidean rotation is nonsingular. In modular flow convention with period 1 in imaginary modular time, the analytical two-point functions are evaluated at Euclidean Killing time βs, 0<s<1. Wick's rule supplies the corresponding connected four-point function.

Choose a nonzero nonnegative v∈C_c^∞(U) in a small interior ball, f_{x,ε} a smooth positive unit-integral spatial mollifier supported at distance O(ε) of each x∈supp v, and h a fixed smooth unit-integral spatial smear supported in a disjoint interior ball. Let F_{x,ε}=f_{x,ε}−h and
```
A_ε = ∫ v(x) (φ(F_{x,ε})² − ω_0(φ(F_{x,ε})²)) dΣ_x.
```
A_ε denotes the L² random variable of the commuting smeared equal-time fields; bounded local tests approximating it will be used below. For fixed ε its state-response derivative d_ε is finite. The order −4 kernel in (2) has near-diagonal coefficient +(8π²)⁻¹ log(R/d) in four Euclidean dimensions (principal Fourier symbol |p|⁻⁴). Off diagonal it is smooth, and integration against fixed h gives bounded contribution. At the equator and uniformly in x on the compact support,
```
d_ε = ∂_μ ω_μ(A_ε)|_{0+}
     = −(1/(8π²))(∫ v dΣ) log(R/ε) + O(1).                            (3)
```
The coefficient follows equally from the flat high-momentum derivative −1/(4|p|³) of the time-zero two-point symbol and radial integration in three dimensions. Curvature, the homogeneous mode, and choice of h affect only the bounded remainder.

Write B_ε for the centered Kubo–Mori covariance of A_ε at ω_0. For fixed ε the modular integral equals its Euclidean form:
```
B_ε = 2 ∫_0^1 ds ∫ v(x)v(y)
          [G_0^⊥(F_{x,ε}, R_{βs}F_{y,ε})]² dΣ_x dΣ_y,                  (4)
```
where R_{βs} is rotation of the second spatial slice in Euclidean static time and each test carries its equatorial integration measure. A fixed h end in any term renders the Green pairing uniformly bounded, since 1/d² is locally integrable in the three-dimensional spatial variables including at τ=0. The local Green function between the mollifiers is bounded by C/[d_Σ(x,y)²+τ²+ε²]+C, with τ=β min(s,1−s), on the coincident charts around s=0 and 1 (positive lapse bounded below). For separations d+τ≫ε this follows from unsmeared elliptic kernel estimates; at d+τ≲ε it follows by scaling the spatial averages of the locally integrable 1/d² singularity. Away from those charts the rotated interior supports have no coincident points and the pairing is bounded. Therefore the relevant radial integral is bounded by
```
C ∫_0^{O(R)} r³ dr/(r²+ε²)² + C = O(1+log(R/ε)).
```
Thus 0≤B_ε≤C(1+log(R/ε)). The imaginary modular integration is essential for this bound: with an equal-time variance instead, the singular squared propagator is integrated in only three relative dimensions and grows more strongly, giving no such relative-entropy witness.

## Entropy lower bound and domain control

For a bounded self-adjoint A∈M with ω_0(A)=0, let P_0(sA) be Araki perturbation pressure. The relative entropy variational inequality and pressure derivatives give
```
S_M(ω_μ||ω_0) ≥ s ω_μ(A)−P_0(sA),
P_0(sA) = s² B_0(A,A)/2 + o(s²).
```
Here B_0 is the connected modular Duhamel covariance, equivalently the Kubo–Mori form. These statements hold in the local type-III algebra, not through a reduced density matrix. If d_A=∂_μω_μ(A)|_{0+} exists, setting s=μ d_A/B_0(A,A) and taking μ↓0 yields a lower bound d_A²/B_0(A,A) on the left side of (1); if the denominator vanishes with nonzero derivative, arbitrary linear test multipliers already give divergence. For positive denominator one first holds A fixed, so no remainder uniformity in its norm or ultraviolet scale is assumed.

To apply this to the unbounded formal A_ε, fix ε first, approximate the spatial integral by finite Riemann sums of the smeared squares, and center. For these quadratic finite sums, the response derivative converges to (3) by smooth covariance differentiation on the fixed smeared family. Their BKM form converges to (4) by L²(ω_0) convergence, or by the Wick integral. To see L² suffices in the modular direction, the KMS Duhamel covariance of a centered self-adjoint test is bounded by its variance (the logarithmic mean ≤ arithmetic mean, or modular-strip Cauchy–Schwarz). For each finite sum, replace the associated real polynomial by bounded spectral truncations and center. The covariance converges by L². The derivative also converges: the finitely many equal-time fields have a nondegenerate Gaussian covariance after removal of linear dependencies, and differentiating its density multiplies it by a constant-plus-quadratic Gaussian score in L². Cauchy–Schwarz controls the truncation error. Hence for each ε a bounded A can realize (3) and (4) to arbitrarily small additive errors in derivative and BKM covariance. This deals with both the perturbation domain and the nonuniformity of the ε→0 limit, without invoking a formal divergent Wick-square as an observable of finite modular norm.

Combining (3) and (4) and varying ε **after** the mass derivative bound gives (1), quantitatively a lower bound increasing at least as log(R/ε) for the family of bounded tests. At every finite spatial resolution held fixed the corresponding test is regular; it is the supremum over measurements in the unchanged static-patch von Neumann algebra that is singular. This method does not require proving that the relative entropy is infinite at each fixed m>0; that stronger statement would need finite perturbation pressure control uniform in ε. In particular it would be inaccurate to infer that state normality or Hilbert–Schmidt Bogoliubov implementability implies a finite information Hessian.

## Primary sources and boundary of claim

* [Kirsten–Garriga, arXiv:gr-qc/9305013](https://arxiv.org/abs/gr-qc/9305013), quantizes the massless de Sitter zero sector and identifies shift-invariant observables. [Fewster–Verch, arXiv:1109.6732](https://arxiv.org/abs/1109.6732), treats the covariant massless current and the rigid-shift degeneracy. Neither supplies the modular susceptibility argument.
* [Bros–Moschella, arXiv:gr-qc/9511019](https://arxiv.org/abs/gr-qc/9511019), develops Euclidean/thermal wedge analyticity; [Gérard–Wrochna, arXiv:1706.08942](https://arxiv.org/abs/1706.08942), provides general Wick-rotated Cauchy covariance/Calderón theory. [Petz, A Variational Expression for the Relative Entropy](https://math.bme.hu/~petz/pdf/35variational.pdf), gives the von Neumann algebra perturbation variational formulation.
* [Casini–Teste–Torroba, arXiv:1611.00016](https://arxiv.org/html/1611.00016), §§2–3 and Appendix A, identifies different-mass theories by common Cauchy data and discusses the spacelike UV threshold Δ<d/2, with separately different null-surface comparisons. The present result is consistent with the marginal four-dimensional mass term; the scoped addition is a direct lower bound on the actual shift-neutral de Sitter static-patch information susceptibility, with the infrared zero sector and local normality treated independently. It is not a new flat-space UV scaling principle.

## Sharpened logarithmic coefficient of the test family

In the bounded-test exhaustion used above one can sharpen (4) before any mass limit is exchanged with a cutoff. Write N(x) for the static lapse relative to central proper time on the equator; it is positive on the compact support of v. The two endpoint strips of the angular parameter s jointly give one complete local Euclidean relative-time integral with measure $dT_E/(\beta N(x))$. In normal local coordinates the time-rotated Green singularity is $(4\pi^2|z|^2)^{-1}$, where z is the four-dimensional relative displacement. In the range $C\epsilon<|z|<c$ the two normalized spatial averages replace this leading term by a relative error $O(\epsilon/|z|)+O(|z|)$, and all h terms are uniformly bounded. The Taylor variations in volume, lapse and v add $O(|z|)$. These relative errors are integrable against the logarithmic radial measure; the inner region $|z|\le C\epsilon$ contributes $O(1)$ by the same smearing scaling as in (4). The sphere area $2\pi^2$ and Wick factor 2 then give
$$
 B_\epsilon=\frac{\log(R/\epsilon)}{4\pi^2\beta}
       \int_U \frac{v(x)^2}{N(x)}d\Sigma_x + O(1).
$$
Combining this with (3), the variational susceptibility witnessed by these bounded approximations has lower coefficient
$$
 \frac{d_\epsilon^2}{B_\epsilon}=
 \frac{\beta}{16\pi^2}\frac{(\int_U v)^2}{\int_U v^2/N}
       \log(R/\epsilon)+O(1).
$$
This is the logarithmic scale of a prescribed test family, not an asymptotic formula for the infinite susceptibility or the entropy at finite mass. The leading coefficient is an asymptotic of the fixed smooth test family before replacing each member by an arbitrarily close bounded approximation, not a bound uniform in the truncation norms. For all compact nonnegative v, the ratio in this formula is bounded by $\int_{\operatorname{supp}v}N$ (weighted Cauchy–Schwarz), and on any fixed interior aperture open set the supremum is $\int N$ over that aperture by approximating v=N with cutoffs from within. The construction still holds epsilon fixed before taking the one-sided mass-coupling derivative.
