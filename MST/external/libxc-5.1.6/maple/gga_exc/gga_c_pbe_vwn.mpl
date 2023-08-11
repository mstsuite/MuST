(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_pbe_vwn_params *params;

  assert(p->params != NULL);
  params = (gga_c_pbe_vwn_params * )(p->params);
*)

$include "lda_c_vwn.mpl"
$ifdef gga_c_pbe_params
params_a_beta  := 0.06672455060314922:
params_a_gamma := (1 - log(2))/Pi^2:
params_a_BB    := 1:
$endif

mgamma := params_a_gamma:
mbeta  := (rs, t) -> params_a_beta:
BB     := params_a_BB:

tp   := (rs, z, xt) -> tt(rs, z, xt):

(* Equation (8) *)
A := (rs, z, t) ->
  mbeta(rs, t)/(mgamma*(exp(-f_vwn(rs, z)/(mgamma*mphi(z)^3)) - 1)):

(* Equation (7) *)
f1 := (rs, z, t) -> t^2 + BB*A(rs, z, t)*t^4:
f2 := (rs, z, t) -> mbeta(rs, t)*f1(rs, z, t)/(mgamma*(1 + A(rs, z, t)*f1(rs, z, t))):

fH := (rs, z, t) -> mgamma*mphi(z)^3*log(1 + f2(rs, z, t)):

f_pbe  := (rs, z, xt, xs0, xs1) ->
  f_vwn(rs, z) + fH(rs, z, tp(rs, z, xt)):

f  := (rs, z, xt, xs0, xs1) -> f_pbe(rs, z, xt, xs0, xs1):

