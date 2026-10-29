(* ::Package:: *)
(* Historical e20 checks for Theorem W3 (all-(N,r) generic rigidity).
   Part 1 — reach: ZSS certificates at N = 24, 28, 32 (meet-in-the-middle exact).
   Part 2 — the W3 mechanism: on the analytic family V = [I_r | eps C] (C symbolic),
   expand every wedge weight to leading order in eps and verify each nonempty
   subset-sum is a NONZERO polynomial in C, for every row — per row class.  Nonzero
   leading polynomial => the subset-sum is not identically zero on the Grassmannian
   => its zero set is measure-zero => generic witnesses are ZSS-free => Lemma W2
   gives rigidity: the analytic-genericity route to ALL (N,r).                    *)

log[s_] := Print["== ", s];

(* ---------- Part 1: reach ---------- *)
drawP[nn_, rr_, wd_] := Module[{V},
   While[True,
    V = RandomInteger[{-wd, wd}, {rr, nn}];
    If[MatrixRank[V] == rr,
     Module[{P = Transpose[V] . Inverse[V . Transpose[V]] . V},
      If[Min[Abs[Flatten[P]]] > 0, Return[P]]]]]];
zssFreeQ[w_] := Module[{h1, h2, s1, s2, set2},
   h1 = Take[w, Ceiling[Length[w]/2]]; h2 = Drop[w, Ceiling[Length[w]/2]];
   If[MemberQ[Rest[Total /@ Subsets[h1]], 0], Return[False]];
   If[MemberQ[Rest[Total /@ Subsets[h2]], 0], Return[False]];
   s1 = Rest[Total /@ Subsets[h1]]; s2 = Rest[Total /@ Subsets[h2]];
   set2 = Association[Thread[(-s2) -> True]];
   ! AnyTrue[s1, KeyExistsQ[set2, #] &]];
certify[nn_, rr_, wd_: 9, tries_: 4] := Module[{P, den, M, sg, ok, bad, t0, w, att = 0},
   While[att < tries, att++;
    t0 = AbsoluteTime[];
    P = drawP[nn, rr, wd];
    den = LCM @@ (Denominator /@ Flatten[P]);
    M = den P; sg = Sign[M]; ok = True; bad = 0;
    Do[w = Table[If[k != i && k != j, M[[i, k]] M[[k, j]] sg[[i, j]], Nothing], {k, nn}];
     If[! zssFreeQ[w], ok = False; bad++], {i, 1, nn - 1}, {j, i + 1, nn}];
    If[ok,
     Print["   Gr(", rr, ",", nn, "): ZSS PASS (witness ", att,
       ")  ** PROVED **   [", Round[AbsoluteTime[] - t0, 0.1], " s]"];
     Return[True],
     Print["   Gr(", rr, ",", nn, "): witness ", att, " inconclusive (", bad,
      " rows)"]]];
   False];
log["Part 1: reach (N = 24, 28, 32)"];
SeedRandom[2026];
certify[24, 12];
certify[28, 14];
certify[32, 16];

(* ---------- Part 2: symbolic nonvanishing on the eps-family ---------- *)
log["Part 2: W3 mechanism — symbolic leading-order subset sums"];
checkFamily[nn_, rr_] := Module[
   {cs, V, G, P, eps, m1, sg, lead, rows, allOK = True, worst = Infinity, cnt = 0},
   cs = Table[Symbol["c" <> ToString[a] <> "x" <> ToString[m]], {a, rr}, {m, nn - rr}];
   V = Join[IdentityMatrix[rr], eps cs, 2];
   G = V . Transpose[V];
   P = Simplify[Transpose[V] . Inverse[G] . V];
   rows = Flatten[Table[{i, j}, {i, nn - 1}, {j, i + 1, nn}], 1];
   Do[Module[{i = rw[[1]], j = rw[[2]], wks, ords, lead0, Tsubs},
     wks = Table[If[k != i && k != j, P[[i, k]] P[[k, j]], Nothing], {k, nn}];
     (* signed weight w_k = P_ik P_kj sign(P_ij): the sign(P_ij) is common to the
        row, irrelevant for zero/nonzero of subset sums: drop it *)
     ords = Exponent[Series[#, {eps, 0, 8}] // Normal, eps, Min] & /@ wks;
     lead0 = MapThread[SeriesCoefficient[#1, {eps, 0, #2}] &, {wks, ords}];
     (* leading order may differ per wedge: a subset-sum vanishes identically only
        if it vanishes at EVERY order; sufficient check: the lowest order present
        in the subset has a nonzero coefficient sum, else escalate order by order *)
     Tsubs = Rest[Subsets[Range[Length[wks]]]];
     Do[Module[{T = t, omin, coef},
       omin = Min[ords[[t]]];
       coef = Total[If[ords[[#]] == omin, lead0[[#]], 0] & /@ t];
       cnt++;
       If[PossibleZeroQ[Simplify[coef]],
        (* escalate: full series sum to higher order *)
        Module[{full = Simplify[Normal[Series[Total[wks[[t]]], {eps, 0, 10}]]]},
         If[PossibleZeroQ[full], allOK = False;
          Print["     !! identically-zero subset sum at row ", rw, " T=", t]]]]],
      {t, Tsubs}]],
    {rw, rows}];
   Print["   (", nn, ",", rr, "): ", cnt, " subset sums checked symbolically — ",
    If[allOK, "ALL nonzero as polynomials in C  ** mechanism VERIFIED **",
     "FAILURES found"]];
   allOK];
checkFamily[5, 2];
checkFamily[6, 3];
log["e20 complete."];
