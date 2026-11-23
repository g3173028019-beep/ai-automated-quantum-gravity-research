# Airy asymptotics at a Lorentzian Regge 4–2 degeneracy

**Abstract.** A reflection-symmetric six-vertex Lorentzian circuit gives a simple zero of the 4–2 factor \(D\) on one light-cone regular branch. Holding the actual boundary lengths fixed, the reduced single-edge Regge action has a nonzero cubic at the coalescence of a flat and a curved stationary point. The localized real contour has a uniform Airy and Airy-derivative approximation. Its two separated Gaussian terms have opposite signed Fresnel phases, and the flat saddle phase flips when \(D\) changes sign. The family is tangent to the real fold discriminant, not a path into an imaginary-saddle regime.

## 1. A Lorentzian 4–2 circuit

Use \((-+++)\), coordinates \((t,x,y,z)\), dimensionless lengths (a physical common length unit \(\ell\) is restored by using \(\hbar_{\rm eff}=8\pi G\hbar/\ell^2\) in the phase). Vertices 0 and 1 are \((0,-a,0,0),(0,a,0,0)\) with \(a>0\) near 1. The shared tetrahedron consists of the four points \(x=0\) with \((t,y,z)\) equal to
\[
 (\sqrt2,\sqrt3,0),\quad(\sqrt2,-\sqrt3,0),\quad
 (-\sqrt2,0,\sqrt3),\quad(-\sqrt2,0,-\sqrt3).
\]
Their barycenter is the origin, in the interior of the tetrahedron and the open edge 01. The resulting two triangulations of the affine bipyramid are the ordinary 4–2 circuit. Both are embedded flat at \(s_f=4a^2\), where \(s=s_{01}\). The boundary apex-to-tetra squared lengths all equal \(d=a^2+1\); tetrahedron same-pair and cross-pair squared lengths are \(12\) and \(-2\), respectively. The full circuit lies on the Lorentzian quadratic sphere centered at 0 at \(a=1\).

The causal checks can be made at \(a=1\) and are open. Internal triangles \(01i\) are spacelike with orthogonal \(x\)- and transverse projected directions of positive squared norm. Internal tetrahedra \(01ij\) have a transverse \(2\times2\) Gram of signature \((-,+)\): determinant \(-24\) for the same pair and \(-3\) for a cross pair, before the positive edge factor. Shared-tetrahedron triangles of three base points have orthogonal directions of squared norms proportional to \(12\) and \(3-8=-5\); the shared tetrahedron itself spans the full Lorentzian transverse affine 3-space. Each triangle containing one apex and two base points has Gram \(s_{ij}(d-s_{ij}/4)\), nonzero for both pair types. A tetrahedron containing one apex and three base points has a spacelike pair-difference direction, orthogonal to the other two vectors. Their Gram determinant is \((a^2-2)(a^2+1)-(a^2+2)^2=-11\) at \(a=1\), nonzero. Every refined 4-simplex has the spacelike edge and three independent transverse vectors, hence signature \((-,+,+,+)\) and nonzero volume. Each coarse 4-simplex likewise has nonzero volume since its apex is out of the shared tetrahedron plane. All edges have nonzero signed length. Thus there is a common fixed real neighborhood in \((s,a)\) of \((4,1)\) with all relevant simplicial Grams nondegenerate and unchanged causal signs. The glued flat normal planes are light-cone regular and the interior wedge light-crossing pattern cannot change in this neighborhood, because none of its projected rays becomes null. Boundary wedge counts are also fixed; the standard Lorentzian boundary term is taken on the same branch. This is not a continuation from Euclidean simplices.

All four refined absolute 4-volumes at the flat edge are \(v_4=a\sqrt2/2\). [Borissova–Dittrich, arXiv:2303.07367](https://arxiv.org/html/2303.07367), Eq. (55), therefore gives \(D^{(4-2)}=(a^2-1)v_4\), with a simple zero. This is not a gauge orbit: all vertices lie on the boundary and its flat realization at fixed boundary lengths sets the internal edge. Off-shell variation along this edge is the sole refined integration direction and the deficit computed below detects it.

## 2. The signed Regge action and its stationary points

Hold the boundary lengths just specified fixed when varying \(s\). Put \(\delta=(s-4a^2)/4\). In a simplex containing edge 01, the projected transverse Gram of its base vertices has diagonal \(1-\delta\), same-pair dot product \(-5-\delta\), and cross-pair dot product \(2-\delta\). Indeed the projection onto the edge measured from its midpoint subtracts \(s/4\) from all entries, since the vertices are equidistant from the apexes. Symmetries act transitively on the four internal hinges \(01i\). Their areas are all
\[
 A(s)=\tfrac14\sqrt{s(4d-s)},\qquad A'(s)=\frac{2d-s}{4\sqrt{s(4d-s)}}.
\]
To find the deficit at \(01i\), project the other transverse vertices perpendicular to vertex \(i\)'s transverse vector. Write \(m=1-\delta\), \(u=-5-\delta\), \(v=2-\delta\). The other vertex of the same pair then has projected squared norm \(m-u^2/m\), either vertex of the other pair has norm \(m-v^2/m\), and their pairwise dots are \(v-uv/m\) (first type of wedge, twice) or \(u-v^2/m\) (second type of wedge, once). All rays are timelike near the central point. The first two wedges join opposite time cones, so each crosses two null rays; the last wedge lies in the same time cone. On the \(\theta^+\) branch of Borissova–Dittrich Eq. (1), their Lorentzian real boost angles are respectively \(+\operatorname{arcosh}r_1\) and \(-\operatorname{arcosh}r_2\); for the first two complex wedges the extra constant in \(\theta^+\) is \(-\pi\) each. Thus those constants cancel the bulk \(2\pi\) in \(iS=\sum_h\sqrt{\mathbb V_h}\epsilon_h\) (their Eqs. (2)–(3)). The **real** deficit factor multiplying each spacelike area in the derivative of \(S\) is
\[
 e(s,a)=\operatorname{arcosh}r_2-2\operatorname{arcosh}r_1,\qquad
 r_1=\frac{\sqrt3(2-\delta)}{\sqrt{(2+\delta)(3-2\delta)}},\quad
 r_2=\frac{9-8\delta}{3-2\delta}.
\]
At \(\delta=0\) these are \(\sqrt2,3\); hence \(e=0\). Direct differentiation gives \(\partial_\delta e=2\sqrt2/3\) there, so \(e_s(s_f,a)=\sqrt2/6\). The Schläfli identity, applied to the full Regge action including its boundary curvature angles, gives exactly \(S_s=4A'(s)e(s,a)\), since no boundary area changes under this variation. In particular at \(a=1,s=4\), \(A'=0,\ A''=-1/16\), and
\[
 S_s=S_{ss}=0,\qquad S_{sss}=8A''e_s=-\frac{\sqrt2}{12}.
\]
In our fixed action sign, the exact local factorization is
\[
 S_s(s,a)=-g(s,a)(s-s_f)(s-s_c),\qquad s_c=2d,\quad
 g=\frac{e(s,a)}{(s-s_f)\sqrt{s(4d-s)}},\quad g(4,1)=\frac{\sqrt2}{24}>0 .
\]
Here \(g\) is analytic via divided differences. There are exactly two nearby stationary points, real and on the same causal branch: the embedded flat point \(s_f\) and the point \(s_c\) with extremal hinge area and nonzero deficit when \(a\ne1\). In particular the second root is not a pure-gauge direction. Writing \(\Delta=s_c-s_f=2(1-a^2)\), the flat Hessian is \(H_f=g(s_f,a)\Delta\), in agreement with the sign opposite to the \(D\) above, and \(H_c=-g(s_c,a)\Delta\). At crossing their **labels**, not the sign ordering of the two sorted real roots, interchange. The exact flat Hessian is $H_f=\sqrt2(1-a^2)/(12a)$. Each coarse 4-volume is $2v_4$, since its determinant is the apex height times the shared tetrahedron determinant. Thus Eq. (78) gives $|H_f|=|D|/(12v_4^2)$, exactly the absolute value obtained here. Its stated flat Fresnel phase $e^{-i\operatorname{sgn}(D)\pi/4}$ also matches.

## 3. Uniform approximation on the real contour

Use the actual edge integration \(ds\) along the real causal branch, and a smooth cutoff equal to 1 on a fixed common neighborhood of both coalescing stationary points, supported inside the verified causal region. Let the amplitude \(m(s,a)\) (including the cutoff) be smooth and nonzero at the coalescence before multiplication by the cutoff. This tests a regular measure, not a measure chosen to divide out the stationary determinant. The resulting uniform amplitude pertains to this pair of saddles in the physical real contour; any far saddle, endpoint, or measure singularity in an uncut global integral is separate and not determined here.

Set \(H=\hbar_{\rm eff}\), let \(s_-<s_+\) be the two sorted saddle values, \(\Phi=[S(s_-)+S(s_+)]/2\), and \(\frac43\zeta^{3/2}=S(s_+)-S(s_-)\). The latter difference is \(\int_{s_-}^{s_+}g(s,a)(s-s_-)(s_+-s)\,ds\), so \(\zeta\) is smooth, real, and
\(\zeta=g(4,1)^{2/3}\Delta^2/4+O(|\Delta|^3)\).
The ordinary real two-saddle cubic normal-form lemma (stationary points and critical values matched, continuous through the nonzero cubic) yields an increasing smooth coordinate \(U\) near the pair with
\[
 S=\Phi-(U^3/3-\zeta U),\qquad
 s_\pm\longleftrightarrow \pm\sqrt\zeta,\qquad
 U_s(4,1)=g(4,1)^{1/3}.
\]
One may apply the analytic fold normal-form theorem in the two-parameter unfolding \((a,p)\) given by adding \(p s\) to the analytic phase and then set \(p=0\). This makes explicit why no singular choice of stationary labels is needed on our tangent path. The push-forward amplitude is \(M(U,a)=m(s(U,a),a)\,ds/dU\). Write
\[
 A_0=\frac{M(\sqrt\zeta)+M(-\sqrt\zeta)}{2},\qquad
 B_0=\frac{M(\sqrt\zeta)-M(-\sqrt\zeta)}{2\sqrt\zeta}
\]
by their smooth divided-difference extensions at zero. Subtracting \(A_0+B_0 U\) in a neighborhood of the two saddles leaves \((U^2-\zeta)\) times a smooth uniform remainder (successive one-variable divided differences); use this factor as a derivative of the cubic phase and integrate by parts. Outside the normal-form neighborhood the local cutoff has no stationary support, and the replacement of the \(A_0+B_0U\) integral by its whole-real-line oscillatory value incurs an \(O(H)\) remainder uniformly. The result is
\[
 \int ds\, m(s,a)e^{iS(s,a)/H}
 =
 2\pi e^{i\Phi/H}
 \left[H^{1/3}A_0\,\mathrm{Ai}(-\zeta/H^{2/3})
      +i H^{2/3}B_0\,\mathrm{Ai}'(-\zeta/H^{2/3})\right]+O(H).
\]
The sign on \(\mathrm{Ai}'\) follows from the negative cubic orientation; the real cycle can equivalently be defined by the ordinary Abel limit. Both terms are needed to recover the two generally unequal stationary amplitudes. At coalescence the leading term instead scales as \(H^{1/3}\) with nonzero coefficient \(2\pi e^{i\Phi/H}\mathrm{Ai}(0)m(4,1)/g(4,1)^{1/3}\).

For \(\zeta/H^{2/3}\gg1\), the Airy oscillatory expansion gives the two ordinary terms \(m(s_j,a)\sqrt{2\pi H/|S_{ss}(s_j,a)|}\exp[iS(s_j,a)/H+i\pi\,\operatorname{sgn} S_{ss}(s_j,a)/4]\). Indeed \(S_{UU}(\pm\sqrt\zeta)=\mp2\sqrt\zeta\), so the two signs and the density Jacobians agree individually. The flat branch reproduces the usual \(|D|^{-1/2}\) singular prefactor only in its separated Gaussian regime; the second real stationary branch cannot be dropped near the crossing on this same real contour. In this family \(\zeta\) remains nonnegative, and the two saddle contributions are on the real contour on both sides. Different complex-cycle prescriptions or a transverse excursion into an unfolding with \(\zeta<0\) carry their own Stokes data; no assertion of universal Lorentzian contour equivalence or refinement invariance is inferred.

## 4. Relation to the Gaussian result

 [Borissova–Dittrich](https://arxiv.org/html/2303.07367) §§II, IV, VII and Appendix A give the causal angle branch, Hessian factor and Gaussian Fresnel convention, and [Dittrich–Kamiński–Steinhaus](https://arxiv.org/abs/1404.5288) treats the six-vertex nonfactorization obstruction. These are inputs; the result above computes the actual surviving Regge vertex and its uniform two-saddle continuation on a verified real Lorentzian family. The statement here concerns the localized contribution of a regular edge density on one causal branch; it neither selects the full path-integral measure nor proves invariance of the uncut move amplitude.
