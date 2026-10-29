# Native mixed quartic interaction and a finite matter-supported sector

29 October 2026, Beijing. Human authors: . Affiliations: .

The starting action is the minimal, massless Einstein--scalar action

\[
 S[g,\phi]=\frac1{2\kappa}\int\sqrt{|g|}R
 -\frac12\int\sqrt{|g|}g^{\mu\nu}\phi_\mu\phi_\nu,
 \qquad \kappa>0.
\]

We retain the historical sign, source and curvature-anchor conventions. The Einstein calculation below uses the integrated Gamma--Gamma representative; no claim about a null-boundary action is needed. The new calculation is a flat, off-shell, physical **mixed quartic contact operation**, including the actual auxiliary elimination. It is distinct from the previously constructed contact induced by a cubic coordinate change. A second calculation constructs an exact finite Einstein--scalar spacetime in the retained Walker sector.

## 1. Actual physical inclusion and native quartic density

Use real transverse coordinates with

\[
 \eta=-2du\,dv+dx^2+dy^2,\qquad
 g_{v\mu}dx^\mu=-du,
 \quad g_{ij}=T_{ij}=\delta_{ij}+C_{ij},
\]

\[
 C=\begin{pmatrix}f&g\\g&-f\end{pmatrix},\quad
 \rho^2=f^2+g^2,\quad \nu=\sqrt{1-\rho^2},\quad
 g_{ui}=B_i,\quad g_{uu}=A.
\]

The physical linear inclusion is

\[
 B_i=\partial_v^{-1}\partial_jC_{ij},\qquad
 A=\partial_v^{-2}\partial_i\partial_jC_{ij}.
\]

In the finite Fourier formulas, \(p=(d,a,p_x,p_y)\) is a covector, \(a\ne0\), and

\[
 B=Cp_\perp/a,\qquad A=p_\perp^TCp_\perp/a^2.
\]

Thus the metric inclusion \(h(p)\) is transverse and trace free for arbitrary external wave symbol \(W(p)=-p^T\eta p\), including \(W(p)\ne0\).

The exact inverse metric gives the scalar density

\[
 \mathcal L_\phi=\nu\left[
 \phi_u\phi_v+\frac A2\phi_v^2
 -\phi_vB^TT^{-1}\nabla_\perp\phi
 -\frac12\nabla_\perp\phi^TT^{-1}\nabla_\perp\phi
 -\frac12B^TT^{-1}B\,\phi_v^2\right].
\]

Since \(C^2=\rho^2I\), its degree \(C^2\phi^2\) term is

\[
\boxed{\mathcal L^{\rm native}_{C^2\phi^2}
 =-\frac{\rho^2}{2}\phi_u\phi_v
 +\phi_v B^TC\nabla_\perp\phi
 -\frac{\rho^2}{4}|\nabla_\perp\phi|^2
 -\frac{|B|^2}{2}\phi_v^2.}
\tag{1}
\]

The \(B\)-terms matter for transverse-dependent waves. Keeping only a transverse metric block would omit them. There is no \(A C\phi^2\) term at this degree: the exact inverse contains \(A\) affinely, and \(\nu A\) has degrees one and three in the physical metric field. The physical fields used in the historical canonical convention are obtained by \((f,g)=\sqrt{2\kappa}(F,G)\); hence (1) is order \(\kappa C_{\rm canonical}^2\phi^2\).

For arbitrary two metric legs \(h_1,h_2\), the same fully polarized native kernel is particularly compact. Define

\[
 D_{12}=\eta h_1\eta h_2\eta+\eta h_2\eta h_1\eta
 -\frac12\operatorname{tr}(\eta h_1\eta h_2)\eta.
\]

Then, for scalar legs 3 and 4,

\[
\boxed{Q^{\rm native}_{12;34}=p_3^TD_{12}p_4.}
\tag{2}
\]

Both (1) and (2) follow directly from the native minimal action, before any nonlinear chart map is applied.

## 2. The off-shell gravitational source is not the on-shell conserved current

Let \(R^{12}_{\mu\nu}\) denote the coefficient of two arbitrary TT physical inclusions in the covariant Ricci tensor. It is computed with the exact quadratic Fourier connection:

\[
 \Gamma^{(1)\rho}_{\mu\nu}(p,h)
 =\frac i2\eta^{\rho\sigma}
 (p_\mu h_{\sigma\nu}+p_\nu h_{\sigma\mu}-p_\sigma h_{\mu\nu}),
\]

\[
 \Gamma^{(2)\rho}_{\mu\nu}(1,2)
 =-(\eta h_1)^\rho{}_{\sigma}\Gamma^{(1)\sigma}_{\mu\nu}(2)
  -(\eta h_2)^\rho{}_{\sigma}\Gamma^{(1)\sigma}_{\mu\nu}(1).
\]

Substitution into \(R=\partial\Gamma-\partial\Gamma+\Gamma\Gamma-\Gamma\Gamma\) supplies the finite expression implemented in `mixed_ricci`. No shell or current-conservation assumption enters it.

The actual Euler variation requires the **raised Einstein tensor density**, rather than just the covariant Ricci tensor. For transverse, trace-free inclusions,

\[
 R^{(1)}_{i,\mu\nu}=-\frac{W_i}{2}h_{i,\mu\nu},\qquad R^{(1)}_i=0.
\]

Consequently

\[
\boxed{\mathcal G^{12}
 =[\sqrt{|g|}G^{\mu\nu}]_{12}
 =\eta R^{12}\eta-\frac\eta2\operatorname{tr}(\eta R^{12})
  +\frac{W_1+W_2}{2}D_{12}.}
\tag{3}
\]

The last term contains the inverse-metric raising corrections and the second-order scalar-curvature correction. It vanishes on the two graviton shells but is generally nonzero off shell. Substituting a conserved on-shell Compton current would delete precisely this term and fail to produce the native off-shell quartic operation.

Use the retained auxiliary coordinates \(y=(t,A,B_x,B_y)\), with \(\delta T=(t/2)I\). Their metric basis matrices \(E_A\) are given explicitly in `aux_basis`. The two-graviton auxiliary source is

\[
\boxed{J^{CC}_{A,12}
 =-\frac1{2\kappa}\mathcal G^{12}:E_A.}
\tag{4}
\]

This formula was checked against a second calculation: the coefficient of \(h_1h_2E_A\) in the Einstein Gamma--Gamma action, with the auxiliary momentum \(-p_1-p_2\). The second calculation includes the linear volume contribution for the trace auxiliary leg; it does not use (3). Four off-shell configurations and all four auxiliary coordinates agree exactly.

One component also admits a closed generic expression. The exact geometric identity

\[
 R_{vv}=-\frac12\operatorname{tr}(T^{-1}T_{vv})
       +\frac14\operatorname{tr}(T^{-1}T_vT^{-1}T_v)
\]

is independent of \(A,B\) and of transverse derivatives. Therefore

\[
 \boxed{J^{CC}_{A,12}
   =\frac{a_1^2+a_2^2+a_1a_2}{4\kappa}
     \operatorname{tr}(C_1C_2),\qquad
 t^{CC}_{12}
   =\frac{a_1^2+a_2^2+a_1a_2}{(a_1+a_2)^2}
     \operatorname{tr}(C_1C_2).}
\tag{4a}
\]

The second formula is the actual nonlinear area/trace backreaction \(y=-H^{-1}J^{CC}\). The direct action calculation additionally verifies (4a) with generic independent longitudinal frequencies and arbitrary two transverse polarizations. Its independence of transverse momenta follows from the displayed exact \(R_{vv}\) identity; it is not inferred by extrapolating those longitudinal tests.

## 3. Complete auxiliary feedback at the mixed quartic order

For scalar covectors \(p=(d,a,p_x,p_y)\), \(q=(e,b,q_x,q_y)\), the actual polarized source is

\[
 J^{\phi\phi}_{pq}=\left(
 -\frac{db+ae}{2},\,-ab,\,a q_x+b p_x,\,a q_y+b p_y\right)^T.
\]

For \(k=(d,a,k_x,k_y)=p_1+p_2\), the retained auxiliary inverse is

\[
 H^{-1}(k)=\kappa\begin{pmatrix}
 0&-4/a^2&0&0\\
 -4/a^2&4d/a^3+2k_\perp^2/a^4&2k_x/a^3&2k_y/a^3\\
 0&2k_x/a^3&2/a^2&0\\
 0&2k_y/a^3&0&2/a^2
 \end{pmatrix}.
\]

Writing the eliminated equation as \(Hy+J^{CC}+J^{\phi\phi}=0\), its mixed fourth-order action is the cross term in

\[
 -\frac12(J^{CC}+J^{\phi\phi})^TH^{-1}(J^{CC}+J^{\phi\phi}).
\]

There is no factor \(1/2\) in the polarized mixed term. The resulting **full mixed quartic contact of the reduced physical action** is

\[
\boxed{Q^{\rm reduced}_{CC\phi\phi}(1,2;3,4)
 =p_3^TD_{12}p_4
  -(J^{CC}_{12})^TH^{-1}(p_1+p_2)J^{\phi\phi}_{34}.}
\tag{5}
\]

This is a finite formula for arbitrary off-shell momenta summing to zero. The scalar sector has no \(C\phi\) auxiliary source around zero scalar background, so there are no omitted mixed pairings of that type. Scalar-exchange and radiative-graviton-exchange diagrams belong to the tree functional obtained **from** the reduced action; they are not part of its contact (5).

Only the already declared nonzero longitudinal pivots \(a_1,a_2,a_1+a_2\) occur in (5). There is no inverse external wave or inverse channel wave in this contact. The historical propagating exchanges retain their physical poles.

One exact off-shell example, in the metric normalization used by the code and with \(\kappa=1\), has

\[
 (W_1,W_2,W_3,W_4)=(-5,-49,-35,-151),
\]

\[
 J^{CC}=(-1289/120,-261/80,153/8,339/20)^T,
\]

\[
 Q^{\rm native}=3161/120,\quad Q^{\rm feedback}=-203/96,
 \quad Q^{\rm reduced}=11629/480.
\]

Every momentum and polarization for this example is stored in `results/native-mixed-quartic.json`.

## 4. Literal curvature-anchored relative-carrier operation

Formula (5) can be evaluated directly on the actual retained carrier, not just on abstract canonical waves. In complex transverse coordinates \(\eta=-2du\,dv+2dz\,dw\), write \(p=(d,a,b,c)\) and

\[
 e_z=(c/a,0,1,0)^T,\qquad e_w=(b/a,0,0,1)^T.
\]

The curvature-normalized inclusion is

\[
\boxed{h_B(p;r_B,t_B)=\kappa a^2r_B\,e_we_w^T
                  +\frac{2}{a^2}t_B\,e_ze_z^T.}
\tag{6}
\]

Use the complex transverse \(\eta\) and the corresponding transformed auxiliary basis in (2)--(5), or transform (6) into \(x,y\) first. In either route, each leg of (6) is TT off shell. Thus (5),(6) define all two-helicity choices of the native mixed quartic operation on the literal \((r_B,t_B,\phi)\) carrier.

The resulting operation is at one specified chart. It supplies the missing native input for a quartic comparison of the matter-assisted Real \(b/c\) representative. It does not prove that comparison: if \(f=L(x+u_2+u_3+\cdots)\), the next defect is

\[
 D_4=S_{4,b}(Lx)-S_{4,c}(x)
     +dS_{3,b}(Lx)[Lu_2]+\tfrac12\langle u_2,Ku_2\rangle.
\]

One must compute its wave-ideal decomposition to solve for \(u_3\), with the actual native quartics included. Replacing \(S_4\) by the coordinate-induced \(\phi^4\) contact alone would not form this equation.

## 5. Exact finite backreaction: real repair of a cubic sign obstruction

This separate spacetime construction uses the retained exact Walker reduction

\[
 g_\theta=-2du\,dv+2dz\,dw
 +\theta_{ww}du^2+2\theta_{vw}du\,dz+\theta_{vv}dz^2,
\]

\[
 E(\theta)=\Box_\eta\theta-\det\operatorname{Hess}_{v,w}\theta,
 \qquad R_{uu}=-E_{ww}/2,\quad R_{uz}=-E_{vw}/2,
 \quad R_{zz}=-E_{vv}/2.
\]

For \(\Psi=\epsilon u\), the full equations require \(E=-\kappa\epsilon^2w^2\), up to affine terms annihilated by the displayed Hessian. A tempting cubic \(\theta=q w^2v\) has determinant \(-4q^2w^2\), so it cannot produce \(+\kappa\epsilon^2w^2\) with real \(q\). This sign problem belongs to that ansatz; it is not a no-go for real backreaction.

The constructive repair is

\[
 H=\alpha w+\beta,\qquad
 \boxed{\theta=
 \kappa\epsilon^2\left(\frac{\alpha w^5}{20}
                       +\frac{\beta w^4}{12}\right)
 +\frac{v^2}{2H},\qquad \Psi=\epsilon u.}
\tag{7}
\]

More generally, \(\theta=f(w)+h(w)v^2/2\) obeys

\[
 \det\operatorname{Hess}\theta
 =h f''+\left(\frac{hh''}{2}-(h')^2\right)v^2.
\]

Choosing \(h=H^{-1}\) cancels the \(v^2\) term exactly; choosing \(f''=\kappa\epsilon^2w^2H\) supplies the required source. Since (7) is independent of \(u,z\), its flat wave operator vanishes.

The metric is explicitly

\[
\boxed{g=-2du\,dv+2dz\,dw
 +\left(\frac{\alpha^2v^2}{H^3}
          +\kappa\epsilon^2w^2H\right)du^2
 -\frac{2\alpha v}{H^2}du\,dz+\frac1Hdz^2.}
\tag{8}
\]

Direct full four-dimensional calculation, independent of the reduced Ricci identity, gives

\[
 \det g=1,\quad R_{\mu\nu}=\kappa\epsilon^2(du)_\mu(du)_\nu,
 \quad (\nabla\Psi)^2=0,\quad\Box_g\Psi=0,\quad R=0.
\]

At \(\epsilon=0\), (8) is a vacuum background. The finite matter-supported correction is exactly

\[
 g(\epsilon)-g(0)=\kappa\epsilon^2w^2H\,du^2.
\]

There is no uncomputed higher power of \(\epsilon\) in this solution. For \(H\ne0\), it is a regular local metric patch, even though \(H=0\) may bound the displayed coordinates. No global extension or singularity classification is asserted.

It is also not a re-expression of an aligned-scalar pp wave when \(\alpha\epsilon\ne0\). The scalar Hessian has nonzero block

\[
 \nabla d\Psi\big|_{u,z}
 =\epsilon\begin{pmatrix}
 -\alpha^2v/H^3&\alpha/(2H^2)\\
 \alpha/(2H^2)&0
 \end{pmatrix},\quad
 \det=-\frac{\epsilon^2\alpha^2}{4H^4}\ne0.
\]

Thus its rank is two. For any pp metric with an aligned scalar \(\Psi=F(U)\) and parallel \(dU\), the Hessian is \(F''(U)dU\otimes dU\), of rank at most one. Rank is invariant under coordinate changes, so the coupled geometry (8) is outside that retained pp family.

The real coordinate section in (8) has signature \((2,2)\). Its complex continuation is a finite analytic Einstein--scalar solution, but it is not the Lorentzian Real completion of the earlier cyclic chart construction. In particular, the standard conjugate-helicity condition on \(z,w\) is not supplied here. These zero-longitudinal/background fields are also outside the \(a\ne0\) Fourier carrier used in (5). The two positive results therefore have deliberately different domains.

The Einstein--Hilbert and aligned-scalar kinetic densities restricted to this Walker ansatz both vanish. The full equations were varied before restriction and independently checked as four-dimensional tensor equations. A multiplier action for the potential equation cannot therefore be called a restriction of the native Einstein--scalar action, or used to infer a BV equivalence from this finite example.

## Sources and attribution

The retained local inputs are `historical Einstein-scalar research/continuation-20260927-round2/matter_transfer.md`, `curvature_anchor_and_descent.md`, `walker_mixed_interactions.md`, and `gravitational_compton_completion.md`; the current continuation's `RESULT.md`, `code/scalar_quartic_and_bv.py`, and `continuation/theory/scalar-quartic-review.md` fix the cubic/induced-quartic baseline. The existing cubic Real \(b/c\) representative is used as context, not silently promoted to quartic/BV completion.

The general integration-out interpretation is established in [Arvanitakis, Hohm, Hull and Lekeu, *Homotopy Transfer and Effective Field Theory I*](https://arxiv.org/abs/2007.07942v2). The heavenly-equation framework is established; [Krasnov and Lipstein, *Pure connection formalism and Plebanski's second heavenly equation*](https://arxiv.org/abs/2412.15779) gives modern context. Neither reference is cited as proving the explicit matter-supported family (7) or the particular off-shell normalization in (3)--(6), which are derived and checked here. The known on-shell Compton amplitude is not counted as a new result.

This continuation was produced by GPT models; the exact finite ansatz repair was proposed in an Astra cross-agent calculation and independently verified here. Earlier work was Fable 5-led; July 2026 onward results were predominantly GPT-produced, as reported by the operator, with DeepSeek and multiple GPT models participating in cross-review and manuscript preparation. Human names and affiliations remain blank. Human expert review has not occurred. Historical dates, this local continuation and public disclosure are distinct. No exhaustive originality claim is made.

> I genuinely do not know how authorship, intellectual credit, or ownership should be attributed for these results.
