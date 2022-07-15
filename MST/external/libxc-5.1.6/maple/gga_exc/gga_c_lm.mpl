(*
 Copyright (C) 2017 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_lm_params *params;

  assert(p->params != NULL);
  params = (gga_c_lm_params * )(p->params);
*)

$define lda_c_vbh_params
$include "lda_c_hl.mpl"

(* Constant after equation (2.25) in Hu-Langreth. It is equal to the
   numerical constant in equation (12) of the Langreth-Mehl paper,
   4.28e-3, once one accounts for the factor 2 from the conversion
   from Rydberg to Hartree *)
lm_J := Pi/(16*(3*Pi^2)^(4/3)):

(* Equation (2.23) in Hu-Langreth *)
lm_d := z -> sqrt(opz_pow_n(z, 5/3) + opz_pow_n(-z, 5/3))/sqrt(2):

(* F parameter, see after eqn (2.25) in Hu-Langreth; this yields
0.26181 which is rounds to 0.262 given by Langreth and Mehl. Hu and
Langreth say the external parameter f = 0.15 for comparison with LM,
but that f = 0.17 is preferable *)
lm_F := 2*sqrt(3)*params_a_lm_f / (2*(3/Pi)^(1/6)):

(* First term in eqn (2.25) in Hu-Langreth *)
lm_t1 := (z, xs0, xs1) ->
  -7/(9*2^(5/3)) * (xs0^2*opz_pow_n(z, 4/3) + xs1^2*opz_pow_n(-z, 4/3)):

(* Second term in eqn (2.25) in Hu-Langreth *)
lm_t2 := (rs, z, xt) ->
  2/lm_d(z) * exp(-lm_F*xt*n_total(rs)^(1/6)) * xt^2:

f := (rs, z, xt, xs0, xs1) ->
  + hl_f(rs, z)
  + lm_J*(lm_t1(z, xs0, xs1) + lm_t2(rs, z, xt))*n_total(rs)^(1/3):
