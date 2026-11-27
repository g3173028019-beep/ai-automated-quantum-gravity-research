# Two reasons for the analytic and angular compatibility hypotheses

**Off-forward values need not determine a forward-supported cut.** At a single mass, the identity kernel on the Bose physical angular Hilbert space has nonnegative coefficients in all partial-wave blocks. Its kernel is supported on coincident directions and the antipodal Bose image. Off-forward pointwise values vanish, but the positive approximate-identity kernel of the manuscript has a nonzero pairing at coincidence: the factor \(q^2\) in the radial density compensates the angular volume and is not a zero of the angular kernel. Discarding the supported component would give a zero real density and a zero dispersive \(B\), even though the forward trace of the identity is infinite. This is a counterexample to replacing the matching of angular test pairings by pointwise agreement off forward. It is not a proposed connected crossing-symmetric amplitude.

**Even crossing and smooth real data alone need not control the shrinking complex discs.** In units where the heavy squared-mass levels are \(S_n=n\), take an integer \(L\ge1\), \(J_n=2n^{L+3}\), and normalized D5 zonal polynomial \(P_J(z)=U_J(z)/(J+1)\). Put
\[
 \rho_n(t)=e^{-n^L}P_{J_n}(1+2t/n),\qquad
 \Delta\mathcal A(s,t)=\sum_{n\ge1}\rho_n(t)
       \left(\frac1{n-s}+\frac1{n-u}\right),\quad u=-s-t.
\]
The scalar spin weights are positive. For \(-A\le t\le0\), \(A<1/10\), all angular arguments are physical and \(|P_J|\le1\); the sum is meromorphic in \(s\) and exactly even under \(s\leftrightarrow u\). On the circles \(|s|=N+1/2\) the remote denominators are bounded below by \(2/5\), and the terms \(n\le |s|/2\) cost \(O(|s|^{-1})\). Hence the full fixed-real-transfer large-circle norm is \(O(|s|^{-1})\). Polynomial Markov derivative estimates in \(J_n\) of every fixed order are suppressed by \(e^{-n^L}\). In particular the real transfer dependence is \(C^\infty\) and the weighted coherent real cut is absolutely integrable.

For any compact band of strictly negative \(a\), let \(t_n=a+32i|a|^3/n^2\). The angle \(\theta(t)=\arccos(1+2t/n)\) obeys \(\theta'(t)=-1/\sqrt{-t(n+t)}\) on the branch through the real angle. On that band \(|\operatorname{Im}\theta(t_n)|\gtrsim n^{-5/2}\) and \(|\sin\theta(t_n)|\lesssim n^{-1/2}\). Using
\[
 P_J(\cos\theta)=\frac{\sin((J+1)\theta)}{(J+1)\sin\theta}
\]
and \(|\sin(x+iy)|^2=\sin^2x+\sinh^2y\) gives
\[
 |\rho_n(t_n)|\gtrsim
 n^{-L-5/2}\exp(cn^{L+1/2}-n^L).
\]
Thus the \(n^{-3}\)-weighted sum of disc suprema diverges. The levels are individually entire in \(t\); it is the coherent **norm** of the continuation that fails.

There are also physical positive D5 tensor cuts with this behavior. Embed the scalar harmonic kernel by the normalized spin-two pair singlet \(\Xi_n=\operatorname{id}/\sqrt5\) in the two physical antipodal screens at the direction \(n\); its kernel is \(P_J(n\cdot m)|\Xi_n\rangle\langle\Xi_m|\), multiplied by the positive level weight. Positivity follows from the spherical addition theorem and the pointwise isometric singlet injection; even \(J\) and symmetric singlet obey the Bose relation. The four-leg projector multiplies its residue by a bounded analytic function of \(t/S_n\) whose value at zero is \(1/25\): in the axial basis the singlet cancels the rotation phases, leaving only the analytic plane factors and channel weights. That factor stays nonzero on the indicated discs at large \(n\). This makes the example a positive physical \(s\)-channel tensor cut, not a full Lorentz-tensor amplitude obeying all crossing equations. In particular it does not satisfy the manuscript's combined analytic-forward-moment and full-crossing hypotheses.
