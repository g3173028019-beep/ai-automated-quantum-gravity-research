# Constrained harmonic action and repair calculations

All displayed quadratic actions use `a0=kappa=1`; the dimensionless coefficients are `lambda_bar=lambda/(kappa a0^2)` and `zeta_bar=zeta/(kappa a0^2)`. Multiply by `a0^2/kappa` to recover the action, use `tau=t/a0`, and use scalar harmonics with `integral_S3 f^2 dvol_gamma=1`. A prime is a tau derivative. The real eight-field norm is `|X|^2+|Z|^2`, represented by a complex four-vector with the real part of its Hermitian inner product.

## Integrated scalar action

Let `b=(2-alpha)/4`, `omega=sqrt(3)`, `Lambda=(6+alpha)/4`, and use the scalar perturbation and ADM variables of RESULT equations (3)–(5). The identities needed for closure are

\[
\Delta(fY)=-(L+3)fY+2dY(\nabla f),
\]

\[
\Delta[dY(\nabla f)]=2LfY-(L-1)dY(\nabla f),
\]

\[
\int|\nabla f|^2=L,\qquad
\int|\nabla\nabla f|^2=L(L-2).
\]

The first two follow by differentiating `Y dot Y=1`, `dY dot dY=gamma` and `nabla_i nabla_j Y=-gamma_ijY`; the last follows by integrating the Bochner identity with `Ric_gamma=2 gamma`. They explicitly exhibit the two target/spatial components that must be mixed in a scalar representation.

The integrated gravitational quadratic Lagrangian is

\[
\mathcal L_g=-3\psi'^2-2L\psi'B-LB^2
+\tfrac32(A+\psi)^2+2L\psi(A+\psi)-L\psi^2
-\tfrac\Lambda2(A+3\psi)^2.
\tag{D1}
\]

The negative conformal kinetic coefficient here is before the lapse and shift constraints. It is not a physical ghost verdict. The four scalar-field amplitudes contribute

\[
\mathcal L_{\rm kin}=\frac b2\{(v'-\omega z)^2+(z'+\omega v)^2
+L[(u'-\omega q-B)^2+(q'+\omega u)^2]
-2\omega BL(z-q)\}
\]

\[
\hspace{5mm}+b\omega(3\psi-A)(z'+\omega v)
+\frac{3b+\alpha}{4}(3\psi-A)^2.
\tag{D2}
\]

The final term includes the canonical clock kinetic density. The `-2 omega BL(z-q)` term is essential: it arises from the background time derivative contracted with `-B grad(delta Phi)`. Omitting it changes the momentum constraint and the physical spectrum.

Spatial gradients and the radial potential give

\[
\mathcal L_{\rm spatial}=-\frac b2\{(L+3)(v^2+z^2)-4L(vu+zq)
+L(L-1)(u^2+q^2)\}
\]

\[
\hspace{5mm}-b(A+\psi)(3v-Lu)-\frac{3b}{4}(A+\psi)^2-\frac{bm}{2}v^2.
\tag{D3}
\]

Varying `A,B` gives RESULT (5), with nonsingular denominator `alpha+6`. Substitution gives

\[
\mathcal L_{\rm red}=\tfrac12x'^TKx'+x'^TDx-\tfrac12x^THx,
\quad G=D-D^T.
\]

The symmetric part of `D` is a total time derivative. The full matrices and unreduced polynomial are saved in `results/quadratic_action.txt` to avoid hiding any stiffness entry behind a claimed reduction.

## Direct covariant equations and exceptional gauge check

With the lapse and shift reconstructed, the original scalar equations are

\[
v''-2\omega z'+(L+m)v-2Lu-6(\psi-A)=0,
\]

\[
u''-2\omega q'+(L-4)u-2v-B'-A-\psi=0,
\]

\[
z''+2\omega v'+Lz-2Lq-\omega(A'-3\psi'-LB)=0,
\]

\[
q''+2\omega u'+(L-4)q-2z-2\omega B=0.
\tag{D4}
\]

The clock equation is `A'-3 psi'-LB=0`. These were obtained directly from the covariant scalar equation, independently of differentiating the action. For example its metric source is

\[
\delta\Box\Phi_0=6(\psi-A)\Phi_0
+(A'-3\psi'-LB)\dot\Phi_0
+(B'+A+\psi)dY(\nabla f)
+2B\,d\dot Y(\nabla f),
\]

with the target rotation included in the last term. Substituting `x''=-K^{-1}(Gx'+Hx)` gives four exact zeros and the exact clock zero in the executable calculation.

For the repaired action put

\[
W=\frac83\bar\lambda b^2L(L-3)u^2
+4\bar\zeta b^2L^2(u+\psi)^2+\tfrac12bLh_*q^2.
\]

The `u` and `q` equations (D4) acquire respectively `(partial_u W)/(bL)` and `(partial_q W)/(bL)`; the radial equations retain their displayed form. All four and the clock identity again reduce to exact zeros with `H -> H+Hess W`.

At `L=3`, `n=(0,-1,0,0,1)` obeys `Kn=Gn=Hn=0`, including the repaired H. Set `psi=0` and compute the four-coordinate dynamics from the principal four-by-four matrices. Substituting those accelerations back into **all five** equations gives exact zeros. This verifies that the omitted psi equation is the gauge-dependent copy, rather than an additional restriction on the physical initial data.

The clock equation is unchanged at linear order: clock-slice displacement has no first-order effect on `E,Q` on this static background, and the background field momentum density is zero. This statement is specific to expansion about the retained background, not a claim that the added operators are clock independent nonlinearly.

## Exact positivity certificates

For `alpha=1/10,m=-207/25`, the leading kinetic minors in the ordering `(v,u,z,q,psi)` are

\[
\frac{19}{40},\quad\frac{361L}{1220},\quad
\frac{6859L}{744200},\quad\frac{130321L^2}{29768000},\quad
\frac{130321L^2(L-3)}{14884000}.
\]

They are positive for every `L>=8`. At the repaired point, the stiffness splits into the real `(v,u,psi)` and imaginary `(z,q)` blocks. Their leading minors are

\[
R_1=\frac{19(1525L+4473)}{61000},
\]

\[
R_2=\frac{361L(613600L^2-456633L-1596861)}{183000000},
\]

\[
R_3=\frac{361L(L-3)}{36600000000}
(259509975L^3-1960056193L^2+3767182800L-534407580),
\]

\[
I_1=\frac{19L}{610},\qquad I_2=\frac{361L^2(L-2)}{24400}.
\tag{D5}
\]

Substitution `L=t+8` makes every coefficient of each expression strictly positive, including the constant term. The executable calculation checks and saves the rational coefficients, proving positivity for all `t>=0`.

At exceptional `L=3`, the four leading repaired stiffness minors are

\[
\frac{21489}{7625},\quad\frac{23064651}{1525000},\quad
\frac{1314685107}{930250000},\quad
\frac{74937051099}{37210000000}.
\]

The corresponding kinetic minors are also saved as positive rationals. This completes the all-scalar-harmonic proof.

For the original model `ell>=4`, the same coefficient method at `L=t+24` proves the two original H blocks positive. For `ell=1,2,3`, exact determinant polynomials in `w=omega_squared` are saved, with Sturm counts. In particular `ell=3` has four positive roots and one negative root; `ell=1` has two real positive roots and a nonreal conjugate pair; `ell=2` has four positive roots and a zero.

The general original determinant is

\[
\det H=\frac{L^3\alpha(L-8)(L-3)(\alpha-2)^4(2L-\alpha-6)}{16(\alpha+6)^2}
\]

\[
\hspace{4mm}\times[2L^2-L\alpha+2Lm-14L-\alpha m+24\alpha-6m-48].
\tag{D6}
\]

This proves the radial-potential obstruction in RESULT without trusting a near-collision numerical eigenvalue.

## Vector reduction

Let `T_i` be a unit divergence-free vector harmonic with rough Laplace eigenvalue `1-n(n+2)` and put `delta Phi/r=exp(i omega tau) dY[(V+iZ)T]`. Use a transverse shift `BT`; the spatial vector metric is set to zero by gauge. Write `d=n(n+2)-3`. The full quadratic Lagrangian before eliminating B is

\[
\mathcal L_V=\frac b2[V'^2+Z'^2+2\omega(VZ'-ZV')-d(V^2+Z^2)]
-bB(V'-2\omega Z)+\tfrac12(b+d/2)B^2.
\]

The gravitational term is `d B^2/4`, obtained from `integral |sym nabla T|^2=d/2`; it vanishes for Killing harmonics at `n=1`. The shift is

\[
B=\frac{2b}{d+2b}(V'-2\omega Z).
\]

Repair contributions are `-2 lambda_bar b^2 d V^2-b h_* Z^2/2`; compression has no transverse first variation. Eliminating B gives RESULT (9), which is positive for all `n>=2`. At `n=1`, V disappears from the reduced action and Z has frequency squared `12+h_*`. Its zero was a spatial Killing gauge direction; the physical relative-rotation mode was never a zero-potential-mass free scalar.

## Tensor correction and high-frequency symbol

For spatial TT `h`, `delta Phi=0`, the original closed sector has unit-normalized quadratic action `1/8[h'^2-(n(n+2)+2b)h^2]`. The new traceless-strain variation is `E^(1)=-b h`; therefore `delta omega_T^2=8 lambda_bar b^2`. Compression and Gram have zero linear TT variation, and the scalar equation sourced by the traceless strain is proportional to the TT divergence and trace, hence vanishes.

The two scalar amplitudes `(u,psi)` have large-L matrices

\[
K_\infty=\begin{pmatrix}c&c\\c&c+2\end{pmatrix},\quad
c=\frac{2(2-\alpha)}{\alpha+6},
\]

\[
H_\infty=K_\infty+\frac{16}3\bar\lambda b^2\begin{pmatrix}1&0\\0&0\end{pmatrix}
+8\bar\zeta b^2\begin{pmatrix}1&1\\1&1\end{pmatrix}.
\]

The generalized eigenvalues give the squared speeds `1.0563364327934` and `4.0944969005399`; the other three scalar speeds tend to one. The vector real branch has speed `sqrt(1+4 lambda_bar b)=sqrt(1.19)` and the imaginary branch tends to one. The formal all-L stability proof includes these high harmonics as a mathematical property of the specified quadratic model, without asserting ultraviolet validity of the EFT.

## Canonical state normalization

The executable computes the covariance from a strictly positive canonical Hamiltonian using principal positive matrix square roots. In `results/construction.json`, each `finite_scalar_band_1_to_3[ell].canonical_covariance` uses the order `(Q_1,...,Q_d,P_1,...,P_d)`, with `d=4` for ell1 and `d=5` for ell2/3. Here

\[
Q=\frac{a_0}{\sqrt\kappa}x,\qquad
P=\frac{a_0}{\sqrt\kappa}(Kx'+Gx/2),\quad [Q_i,P_j]=i\delta_{ij}.
\]

Thus the physical dimensionless coordinate covariance is `(kappa/a0^2) C_QQ`. There is no extra sphere-volume factor because harmonics have unit L2 norm. Centered finite-band quantum rms amplitudes shrink as `sqrt(kappa)/a0`; no all-time nonlinear Gaussian-containment claim follows from that scaling.

The evolution CSV uses the **classical dimensionless** coordinates `x,x'`, then the reconstructed harmonic metric amplitudes `A,B` and the dimensionless energy E. It was generated by exact matrix exponentials, not an explicit-Euler method. The pointwise metric additionally multiplies the chosen normalized harmonic and its gradient. The energy bound refers to harmonic coordinate norm, not a gauge-invariant pointwise norm on all metric components.
