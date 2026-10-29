# Four positive-energy local clocks for the weighted Maxwell echo

29 October 2026, Beijing. Units: the retained dimensionless model units, with hbar=1. The construction enlarges the completed reduced Maxwell/Bianchi Hamiltonian; it does not derive a covariant apparatus stress tensor.

## 1. Retained system and the particular obstruction being repaired

Let H be the self-adjoint, completed finite-code Hamiltonian of the retained weighted echo. Its trusted-region expression is

\[
H=-Ap_q^2+B_{\rm sh}q^{-2}P_\beta^2
  +\delta(F_A+F_B)q^{-2/3}S(\beta)+\zeta q^{-2}.
\]

Retain A=1/2, B_sh=2/9, delta=.03, Cmax=.06, |zeta|<=10, the Q-dependent smooth exterior completion and its global matter bound

\[
\|V_m\|\le V_*=12C_{\max}Q^{-2/3},\qquad
V_*/E_0=48k\sqrt\nu.
\tag{1}
\]

Here nu=AT/Q², E0=1/(4nu T), k=2 sqrt(pi/b), and b=.005512414487283133... is the retained echo coefficient. Two simultaneous local X pulses centered at t_j=x_jT, with

\[
(x_1,x_2)=(.1826174653\ldots,.6744044056\ldots),
\tag{2}
\]

close both weighted conditional displacements while retaining |Phi|=pi and the physical nonnegative Maxwell coefficient delta[1+y(F_A+F_B-1)]. The common identity force is retained. The target is a Bell orbit from |++>, up to known party-local phases.

The previous driven construction used a joint compression P H_c P throughout each pulse to keep the background-H support negative. Here that control is removed. Only the initial total-energy resource below is spectrally specified.

A single free quadratic clock that must remain narrow at both events has a genuine dispersion/energy tradeoff over their macroscopic separation. Applying the old narrow pulse width to such a clock is not justified. The repair is four independent one-shot clocks: one clock for each party and each event. Each clock focuses at its own event and is permitted to be broad at other times.

## 2. Explicit time-independent, party-local apparatus

For r=(party j, event l), introduce an internal clock coordinate Q_r on the whole real line and its canonical P_r. These are apparatus coordinates, distinct from the radial geometry q. Choose M>0 and a positive central clock momentum p0, with v=p0/M. Define

\[
f(s)=\begin{cases}1+\cos(2\pi s),&|s|\le1/2,\\0,&|s|>1/2,\end{cases}
\quad \int f(s)ds=1,\quad 0\le f\le2,\quad \operatorname{Lip}(f)=2\pi,
\]

\[
g(Q)=\frac{\pi}{2\tau}f\!\left(\frac{Q}{v\tau}\right),
\quad g_{\max}=\frac\pi\tau,
\quad h_r=\frac{P_r^2}{2M}+g(Q_r)X_j+g_{\max}I.
\tag{3}
\]

The whole apparatus term is positive as a quadratic form:

\[
h_r\ge\frac{P_r^2}{2M}+g_{\max}-g(Q_r)\ge0.
\tag{4}
\]

The offsets have been included in the total-energy ledger, rather than silently choosing a negative apparatus energy zero. In particular the energy resource is P²/(2M), never an unbounded-sign linear momentum clock. Equation (3) by itself is a reduced kinetic clock, with no physical rest mass assigned. The rest-energy-inclusive point in section 7 adds a positive Mc_light² to each h_r and charges it explicitly. The total autonomous Hamiltonian is

\[
\boxed{\mathcal H=H+\sum_{r=1}^{4}h_r.}
\tag{5}
\]

Its controls depend only on each party's own clock and X_j, contain no geometry operators and no A-B control term, and have no externally prescribed time argument. Boundedness of g makes (5) a bounded perturbation of the self-adjoint tensor sum of H and four positive free kinetic operators. Two clocks at a party interact only through that party's qubit. The independent initial clock packets encode the agreed schedule.

## 3. Initial focused packets and their dispersion cost

Let each central kinetic energy be E_c=p0²/(2M)=alpha E0. Choose momentum standard deviation sigma_p=epsilon p0. At its event t_l, a clock's freely evolved packet is the minimum-uncertainty Gaussian centered at Q=0 with mean p0 and

\[
\sigma_*\equiv\sigma_Q(t_l)=\frac1{2\epsilon p_0}.
\]

An explicit initial momentum wavefunction is

\[
\widetilde\chi_l(p,0)=
 (2\pi\sigma_p^2)^{-1/4}
 \exp[-(p-p_0)^2/(4\sigma_p^2)]\,
 \exp[+ip^2t_l/(2M)].
\tag{6}
\]

The plus sign in the chirp is essential: free evolution cancels it at t_l. Both parties use independent copies. Under free propagation,

\[
\bar Q(t)=v(t-t_l),\quad
\sigma_Q(t)^2=\sigma_*^2+\frac{\sigma_p^2}{M^2}(t-t_l)^2,
\quad \operatorname{Cov}(Q,P)_{t=0}=-t_l\sigma_p^2/M.
\tag{7}
\]

Thus the packet is broad and correlated initially, rather than a narrow packet propagating without dispersion. Setting s=(t-t_l)/tau gives the exact free position distribution

\[
\frac{Q(t)}{v\tau}\sim N(s,d(s)^2),\quad
d(s)^2=a^2+\epsilon^2s^2,\quad
a=\frac{\sigma_*}{v\tau}=\frac1{4\epsilon\alpha c},\quad c=E_0\tau.
\tag{8}
\]

The focus timing spread is a tau, and the spread at a window edge is tau sqrt(a²+epsilon²/4). The energy price of focusing is explicit:

\[
\langle K_r\rangle=\alpha E_0(1+\epsilon^2),\quad
\operatorname{sd}(K_r)=\alpha E_0\sqrt{4\epsilon^2+2\epsilon^4}.
\tag{9}
\]

The construction costs highly correlated initial position-momentum preparation and macroscopic position coherence. It makes no claim that this preparation is simple or laboratory ready.

## 4. Full clock/control error, including feedback and reflection

Compare (5) on the product seed to the time-dependent reference in which each clock propagates freely and g(Q_r) is replaced by g(v(t-t_l)). The matter reference is exactly the bare local finite-pulse protocol. Its clock state remains a product, so Duhamel's formula yields a vector-norm error bounded by

\[
\delta_{\rm clock}\le\sum_{r=1}^4\frac\pi2\int ds\ e(s),\quad
e(s)=\left\{\mathbb E\left[
f(s+d(s)Z)-f(s)\right]^2\right\}^{1/2},\quad Z\sim N(0,1).
\tag{10}
\]

This compares the full interacting autonomous evolution with the free-clock reference. It already bounds quantum clock feedback, phase/area errors, wavepacket distortion and reflected amplitude; treating the exact clock as free has not been assumed. No separate heuristic backreaction error is added or omitted.

For |s|<=1, Lipschitz continuity and sqrt(a²+epsilon²s²)<=a+epsilon|s| give

\[
\int_{-1}^1e(s)ds\le2\pi(2a+\epsilon).
\tag{11}
\]

For a<=epsilon, |s|>=1, the Gaussian density bound on the support [-1/2,1/2] gives

\[
\Pr(|s+dZ|\le1/2)
 \le\frac{1}{\sqrt{2\pi}\epsilon|s|}
      e^{-1/(16\epsilon^2)},
\quad
e(s)\le2(2\pi)^{-1/4}\epsilon^{-1/2}|s|^{-1/2}
e^{-1/(32\epsilon^2)}.
\]

Each actual finite integration range is contained in |s|<=Smax=T/tau. Therefore

\[
\boxed{\delta_{\rm clock}\le
4\pi^2(2a+\epsilon)+
16\pi(2\pi)^{-1/4}\epsilon^{-1/2}\sqrt{T/\tau}
e^{-1/(32\epsilon^2)}.}
\tag{12}
\]

The outer tail is retained analytically, even when machine exponentiation underflows. This is a finite-T bound; it is not a uniform infinite-time scattering theorem.

There is also a useful conservative negative-momentum diagnostic. For free clocks the probability that any P_r<0 is at most p_-=4 exp[-1/(2epsilon²)]. For the spectrally prepared state below, the vector distance from the product seed is at most sqrt(2) D_src. At any time up to T, the exact negative-momentum probability is consequently bounded by

\[
\Pr_{\rm actual}(\exists r:P_r<0)
\le\left(\sqrt{p_-}+\sqrt2D_{\rm src}+\delta_{\rm clock}\right)^2.
\tag{13}
\]

This includes initial projection and interaction-induced reflected/negative-momentum amplitude. It does not claim an asymptotic S-matrix reflection coefficient. No mean of an unbounded momentum or kinetic-energy observable is inferred from trace-distance closeness. At the chosen point gmax/Ec=6.3e-7, while the incident relative momentum spread is .001; the rigorous bound remains (12)-(13), rather than a frozen-X classical barrier approximation.

## 5. Finite pulses and the retained nonlinear dynamics

The deterministic pulse (3) has support length tau and angle pi/2 per party, hence implements -iX_j in the absence of drift. Its even profile places the half-angle at the specified switching center. Passing into its local toggling frame leaves the geometry Hamiltonian unchanged. Within each of the two windows, the difference between the rotated matter potential and the instantaneous-echo potential is bounded by V_*. Therefore the shape of the raised-cosine pulse does not change the old bounded-drift comparison:

\[
\delta_{\rm finite}\le2\tau V_*=96ck\sqrt\nu.
\tag{14}
\]

The retained continuum comparison, including the common identity force, gives

\[
\eta_{\rm dyn}\le K_\beta\nu^{1/4}+K_r\sqrt\nu+640\nu,
\quad K_\beta<344.301,\quad K_r<2173469.
\tag{15}
\]

These are the fresh round-02 nonlinear constants, not the older unrefocused propagation coefficients. The accepted specific radial Gaussian and compact shape source are reused. At |Phi|=pi,

\[
Q=\left[\frac{\pi A^2}{b\delta^2\nu^3}\right]^{3/8},\quad
T=\nu Q^2/A,\quad w=\nu^{3/8}.
\tag{16}
\]

The true autonomous comparison uses (12) after (14)-(15); the ideal weighted-loop phase is not simply asserted to survive a quantum clock.

## 6. One initial joint resource, then exact dust positivity forever

Let Psi_s be the product of the retained geometric/shape seed, |++> and the four clocks (6). The retained seed graph bound is

\[
\|(H+E_0)\psi_s\|\le R_0E_0\sqrt\nu,
\quad R_0=10+50k<2397.284.
\]

Charge all four positive offsets and define

\[
\beta=1-4\alpha(1+\epsilon^2)-4\pi/c-\rho_{\rm rest}>0,
\quad
r=R_0\sqrt\nu+2\alpha\sqrt{4\epsilon^2+2\epsilon^4}+4\pi/c.
\tag{17}
\]

Here rho_rest=0 for the pure reduced kinetic model, or rho_rest=4Mc_light²/E0 when each rest energy is included. Rest energies are constants and add no graph variance. The first clock contribution to r is the standard deviation of four independent kinetic energies, not the sum of their means. The final term bounds the initial bounded interactions globally, without relying on absent initial packet tails. The source obeys

\[
\|(\mathcal H+\beta E_0)\Psi_s\|\le rE_0.
\]

Specify the prepared initial resource by

\[
P_{\rm aut}=\mathbf1_{(-\infty,-\beta E_0/2]}(\mathcal H),\quad
\Psi_0=P_{\rm aut}\Psi_s/\|P_{\rm aut}\Psi_s\|.
\tag{18}
\]

By the spectral theorem, the discarded amplitude is at most 2r/beta. For pure states its square is the source failure probability, while the trace distance from Psi_s is exactly the discarded amplitude:

\[
D_{\rm src}\le2r/\beta,\quad
p_{\rm fail}\le(2r/\beta)^2.
\tag{19}
\]

There is no additional normalization factor in this trace-distance formula. Vector-norm comparisons use the separate sqrt(2) factor in (13).

Because the total Hamiltonian is time independent, it commutes with P_aut and preserves its support exactly:

\[
\boxed{-\mathcal H\ge\beta E_0/2>0\quad\text{on the actual evolving state}.}
\tag{20}
\]

Thus dust energy is strictly positive for the enlarged constrained model even during the actual clock interactions. Neither negative expectation alone nor proximity to a negative-energy seed is used to infer (20). No repeated postselection or spectral control occurs during propagation. The same graph bound gives

\[
\|\mathcal H\Psi_0\|\le
\frac{(\beta+r)E_0}{\sqrt{1-(2r/\beta)^2}},
\tag{21}
\]

which is conserved exactly.

Preparation (18) is still a joint geometry/matter/apparatus spectral resource. It is a scientific assumption, not a constructed local preparation instrument or a gate factory acting on unknown qubits. Its entangling contamination is at most D_src and is included below. If one uses only the product seed, the same apparatus and gate comparison work without this source charge, but exact dust support is replaced by a negative-energy failure probability bounded by (19).

## 7. A finite entangling point and resource ledger

Take nu=1e-28, alpha=.05, epsilon=.001, c=1e8, M=1. Then

\[
Q\simeq2.82\times10^{33},\quad T\simeq1.59\times10^{39},
\quad E_0\simeq1.58\times10^{-12},\quad
\tau\simeq6.35\times10^{19},\quad \tau/T=4\times10^{-20}.
\]

All four apparatus terms are nonnegative. Their total mean kinetic cost is .2000002 E0; the four offsets add less than 1.26e-7 E0. Their combined free kinetic standard deviation is about .0002 E0. The prepared total-H RMS is below .801 E0. The exact dust gap is above .399 E0 and the source failure bound is below 2.51e-7.

Each clock has p0 approximately 3.97e-7 and sigma_p approximately 3.97e-10. The focused position width is approximately 1.26e9. The focus timing width is 5e-5 tau and the window-edge width is below 5.03e-4 tau. At t=0 the two event families instead have position widths approximately 1.15e29 and 4.25e29; their centers are approximately -1.15e32 and -4.25e32. These large initial widths and the exact chirp (6) are load-bearing resources, not omitted preparation errors.

| Trace-distance contribution | Conservative evaluated bound |
|---|---:|
| Genuine shape/shear | 0.000034431 |
| Radial nonlinearity/cutoff | 0.000000022 |
| Ordering | 6.4e-26 |
| Initial total spectral resource | 0.000500315 |
| Finite local pulse duration | 0.004583585 |
| Clock timing, dispersion and interaction backreaction | 0.043426260 |
| Total | **less than 0.04855** |

The outer clock-tail term has natural logarithm below -31220 and is included as a strictly positive analytic remainder. Thus, after tracing out both geometry and all four clocks,

\[
\boxed{D(\rho_{AB},\rho_{\rm Bell\ orbit})<.04855,
\qquad N(\rho_{AB})>.451.}
\tag{22}
\]

The Bell-orbit partial-transpose witness has spectral range one, so its expectation shifts by at most the trace distance and yields N>=.5-D. This certifies entanglement in the reduced original matter qubits, without retaining the clocks in the party bipartition. The bound (13) on any negative clock momentum is below .00195 at any fixed time up to T. The Gaussian integral in (10), evaluated independently at Hermite orders 64 and 128, gives a much smaller diagnostic clock error around .00645; the theorem uses the conservative analytic .04343 instead.

Relative to the former nu=1e-16 compressed-control point, this displayed autonomous construction increases Q by about 3.16e13 and T by 1e15. The old compressed-pulse error is removed, while four positive clock resources, one initial total-H filter and a new clock error are charged. It is an existence construction with explicitly large resources, not a claim that autonomous local control is resource free.

### An explicit point that also pays all four rest energies

To remove the ambiguity between a positive reduced kinetic Hamiltonian and a physical mass-energy ledger, choose a second finite point with c_light=1,

\[
\nu=10^{-42},\quad\alpha=10^{-6},\quad
v=.01c_{\rm light},\quad\epsilon=.001,\quad c=10^{14},
\quad M=\frac{2\alpha E_0}{v^2}.
\]

Replace every h_r by Mc_light²+h_r. Then all four rest energies cost

\[
\rho_{\rm rest}=\frac{8\alpha c_{\rm light}^2}{v^2}=.08,
\]

while the four mean kinetic energies cost about 4.000004e-6 E0 and the positive offsets cost less than 1.26e-13 E0. This is still the explicit **Newtonian quadratic-dispersion model**, now with a complete rest-energy term. No relativistic-dispersion approximation over the entire propagation interval is asserted: a small v/c alone does not control its accumulated phase over this enormous T.

At this point,

\[
Q\simeq1.58\times10^{49},\quad T\simeq5.02\times10^{56},\quad
E_0\simeq4.98\times10^{-16},\quad\tau\simeq2.01\times10^{29},
\quad M\simeq9.96\times10^{-18},\quad p_0\simeq9.96\times10^{-20}.
\]

The dimensionless focus width is a=2.5e-6, giving the conservative clock term below .039675810. Finite pulses cost less than .000458359, shape/shear less than 1.09e-8, and the initial spectral source less than 8.70e-9. The outer-tail logarithm remains below -31211. Thus this second finite point achieves

\[
\boxed{D(\rho_{AB},\rho_{\rm Bell\ orbit})<.04014,
\quad N(\rho_{AB})>.459,\quad
-\mathcal H>.459E_0\ \text{on the actual support}.}
\tag{23}
\]

The actual evaluated dust gap is approximately .459998 E0, source failure probability is below 7.57e-17, and total-H RMS is below .921 E0. Equation (13) bounds any negative clock momentum by .001575. The initial packet centers are about -9.17e53 and -3.39e54 and their standard deviations about 9.17e50 and 3.39e51. The focus position width is about 5.02e21; sigma_time at focus is 2.5e-6 tau. The smaller energy price has been paid with longer pulses, still smaller nu, and much larger geometric and initial clock resources. It has not been obtained by discarding rest energy.

## 8. Scope

This is an autonomous, self-adjoint enlarged reduced quantum model with positive apparatus energy, exact positive-dust support on a specified initial resource, and a controlled reduced-matter entangling output. The second point also explicitly pays rest energies. It replaces the nonlocal dynamical pulse compression. It does not replace the initial spectral preparation by a local instrument, establish relativistic clock dispersion or a fully covariant clock stress tensor, supply full photon QFT, or establish a laboratory realization. The full constrained model has changed by adding apparatus; positivity refers to -H_aut, not to the original background H separately. The source family, completion and resources change with nu.

No new global-priority or originality judgment follows. The accepted echo mechanism and constants remain unchanged; the substantive advance is the four independently focused, positive-energy apparatus construction with a full interacting error and energy ledger.
