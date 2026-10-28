(* ::Package:: *)
(* Historical d11 verification (q=2 profile proof + F(rho) + E_K identity). Self-contained.
   The general sector-variance formula is prior work; see the manuscript attribution.
   T1: exact sector norms <W_{2k}^2> for z = Pf(J), n = 2..5, vs (2n)! C(2m,m)/4^m, m=n-k.
   T2: EGF consistency: x=1 -> <z^4> = (2n+1)((2n-1)!!)^2 ; x=0 -> <z^2>^2.
   T3: E_K det(M + cK) = Sum_r c^{2r} (2r-1)!! Sum_{|S|=2n-2r} det M_S  (symbolic, n=1,2).
   T4: F(rho) = E[det J det(rho J + Sqrt[1-rho^2] K)] = (2n)! Sum_m C(2m,m) 4^-m rho^{2(n-m)}:
       symbolic n=2; Monte Carlo n=3,4. *)

log[s_] := Print["== ", s];

(* ---------- matchings with Pfaffian signs ---------- *)
matchings[verts_List] := Module[{v, rest, res = {}},
   If[verts === {}, Return[{{{}, 1}}]];
   v = First[verts]; rest = Rest[verts];
   Do[Module[{u = rest[[i]], rem, sub},
     rem = Delete[rest, i];
     sub = matchings[rem];
     Do[AppendTo[res, {Prepend[s[[1]], {v, u}], (-1)^(i + 1) s[[2]]}], {s, sub}]],
    {i, Length[rest]}];
   res];

testPf[n_] := Module[{m = matchings[Range[2 n]], J, pf, dt},
   J = Table[0., {2 n}, {2 n}];
   Do[Module[{r = RandomReal[{-1, 1}]}, J[[i, j]] = r; J[[j, i]] = -r],
     {i, 2 n}, {j, i + 1, 2 n}];
   pf = Total[#[[2]] Product[J[[e[[1]], e[[2]]]], {e, #[[1]]}] & /@ m];
   dt = Det[J];
   Abs[pf^2 - dt]];
log["Pf engine sanity (should be ~1e-15):"];
Print["   n=2: ", testPf[2], "   n=3: ", testPf[3]];

(* ---------- T1 ---------- *)
log["T1: exact sector norms, q=2"];
sectorNorms[n_] := Module[{m, idx, c0 = 0, masks, signs, L, recs, grouped, monrecs, monAssoc, norms},
   m = matchings[Range[2 n]];
   idx = <||>;
   Do[c0++; idx[{i, j}] = c0, {i, 2 n}, {j, i + 1, 2 n}];
   masks = Table[Total[2^(idx[#] - 1) & /@ mm[[1]]], {mm, m}];
   signs = m[[All, 2]];
   L = Length[m];
   (* step 1: records {Smask, Umask} -> coeff, grouped *)
   recs = Reap[
      Do[Sow[{BitAnd[masks[[i]], masks[[j]]], BitXor[masks[[i]], masks[[j]]],
          signs[[i]] signs[[j]]}], {i, L}, {j, L}]][[2, 1]];
   grouped = Normal@GroupBy[recs, (#[[1 ;; 2]] &) -> (#[[3]] &), Total];
   (* step 2: expand T over subsets of S bits *)
   monrecs = Reap[
      Do[Module[{S = g[[1, 1]], U = g[[1, 2]], cc = g[[2]], bits, nb},
        If[cc =!= 0,
         bits = Flatten@Position[Reverse@IntegerDigits[S, 2, 64], 1];
         nb = Length[bits];
         Do[Module[{sel, T},
            sel = If[nb == 0, {},
              Flatten@Position[Reverse@IntegerDigits[sub0, 2, Max[nb, 1]], 1]];
            T = If[sel === {}, 0, Total[2^(bits[[sel]] - 1)]];
            Sow[{T, U, cc}]],
           {sub0, 0, 2^nb - 1}]]],
        {g, grouped}]][[2, 1]];
   monAssoc = Normal@GroupBy[monrecs, (#[[1 ;; 2]] &) -> (#[[3]] &), Total];
   (* step 3: norms by chaos order *)
   norms = <||>;
   Do[Module[{T = g[[1, 1]], U = g[[1, 2]], cc = g[[2]], ord, w},
      If[cc =!= 0,
       ord = 2 DigitCount[T, 2, 1] + DigitCount[U, 2, 1];
       w = cc^2 2^DigitCount[T, 2, 1];
       If[KeyExistsQ[norms, ord], norms[ord] += w, norms[ord] = w]]],
     {g, monAssoc}];
   norms];

Do[Module[{t, nr, ok},
   {t, nr} = AbsoluteTiming[sectorNorms[nn]];
   ok = And @@ Table[
      Lookup[nr, 2 k, 0] === (2 nn)! Binomial[2 (nn - k), nn - k]/4^(nn - k), {k, 0, nn}];
   Print["   n=", nn, " (", Round[t, 0.1], "s)  computed: ",
     Normal[KeySort[nr]], "   MATCH vs (2n)!C(2m,m)/4^m: ", ok]],
  {nn, {2, 3, 4, 5}}];

(* ---------- T2 ---------- *)
log["T2: EGF consistency (symbolic)"];
gf[x_, y_] := 1/((1 - x^2 y^2) Sqrt[1 - y^2]);
Do[Module[{z4, z4ref, z2sq},
   z4 = (2 nn)! SeriesCoefficient[gf[1, y], {y, 0, 2 nn}];
   z4ref = (2 nn + 1) ((2 nn - 1)!!)^2;
   z2sq = (2 nn)! SeriesCoefficient[gf[0, y], {y, 0, 2 nn}];
   Print["   n=", nn, ": sum = ", z4, " vs <z^4> = ", z4ref, " (", z4 === z4ref,
     ");  k=0 term = ", z2sq, " vs <z^2>^2 = ", ((2 nn - 1)!!)^2, " (",
     z2sq === ((2 nn - 1)!!)^2, ")"]],
  {nn, {2, 3, 4, 5, 8}}];

(* ---------- helpers for T3/T4 ---------- *)
gaussE[poly_, vars_] := Module[{ex = Expand[poly]},
   If[Head[ex] === Plus, Total[gaussE[#, vars] & /@ List @@ ex],
    Module[{res = 1, d},
      Do[d = Exponent[ex, v];
        If[OddQ[d], res = 0];
        If[EvenQ[d] && d > 0, res = res (d - 1)!!], {v, vars}];
      res*(ex /. Thread[vars -> 1])]]];
antisym[sym_, n_] := Module[{M = Table[0, {2 n}, {2 n}]},
   Do[M[[i, j]] = sym[i, j]; M[[j, i]] = -sym[i, j], {i, 2 n}, {j, i + 1, 2 n}]; M];

(* ---------- T3 ---------- *)
log["T3: E_K det(M+cK) identity (symbolic)"];
Do[Module[{M, K, kvars, lhs, rhs},
   M = antisym[mm, nn]; K = antisym[kk, nn];
   kvars = Flatten[Table[kk[i, j], {i, 2 nn}, {j, i + 1, 2 nn}]];
   lhs = gaussE[Det[M + c K], kvars];
   rhs = Total[Map[Module[{S = #, r},
        r = (2 nn - Length[S])/2;
        c^(2 r) (2 r - 1)!! If[Length[S] == 0, 1, Det[M[[S, S]]]]] &,
      Select[Subsets[Range[2 nn]], EvenQ[Length[#]] &]]];
   Print["   n=", nn, ": identity holds: ", Simplify[Expand[lhs - rhs]] === 0]],
  {nn, {1, 2}}];

(* ---------- T4 ---------- *)
log["T4: F(rho)"];
Fpred[n_, r_] := (2 n)! Sum[Binomial[2 m, m]/4^m r^(2 (n - m)), {m, 0, n}];
Module[{J, K, jv, kv, lhs, rhs},
  J = antisym[jj, 2]; K = antisym[kk, 2];
  jv = Flatten[Table[jj[i, j], {i, 4}, {j, i + 1, 4}]];
  kv = Flatten[Table[kk[i, j], {i, 4}, {j, i + 1, 4}]];
  lhs = gaussE[gaussE[Expand[Det[J] Det[rho J + cc0 K]], kv], jv] /. cc0 -> Sqrt[1 - rho^2];
  rhs = Fpred[2, rho];
  Print["   n=2 symbolic: F = ", Expand[lhs], "   pred = ", Expand[rhs],
    "   MATCH: ", Simplify[lhs - rhs] === 0]];
SeedRandom[42];
mcF[n_, r_, m_] := Module[{c = Sqrt[1 - r^2], acc = 0.},
   Do[Module[{A = Table[0., {2 n}, {2 n}], B = Table[0., {2 n}, {2 n}], x},
     Do[x = RandomVariate[NormalDistribution[]];
      A[[i, j]] = x; A[[j, i]] = -x;
      x = RandomVariate[NormalDistribution[]];
      B[[i, j]] = x; B[[j, i]] = -x, {i, 2 n}, {j, i + 1, 2 n}];
     acc += Det[A] Det[r A + c B]], {m}];
   acc/m];
Do[Module[{est, pred},
   est = mcF[nn, 0.6, 400000];
   pred = N[Fpred[nn, 6/10]];
   Print["   n=", nn, " MC(rho=0.6, 4e5 samples): ", est, "   pred: ", pred,
     "   rel.dev: ", Abs[est - pred]/Abs[pred]]],
  {nn, {3, 4}}];

log["d11 complete."];
