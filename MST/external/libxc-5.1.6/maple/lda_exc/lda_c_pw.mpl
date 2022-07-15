(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_pw_params *params;

  assert(p->params != NULL);
  params = (lda_c_pw_params * )(p->params);
*)

$ifdef lda_c_pw_params
params_a_pp     := [1,  1,  1]:
params_a_a      := [0.031091, 0.015545, 0.016887]:
params_a_alpha1 := [0.21370,  0.20548,  0.11125]:
params_a_beta1  := [7.5957, 14.1189, 10.357]:
params_a_beta2  := [3.5876, 6.1977, 3.6231]:
params_a_beta3  := [1.6382, 3.3662,  0.88026]:
params_a_beta4  := [0.49294, 0.62517, 0.49671]:
params_a_fz20   := 1.709921:
$endif

$ifdef lda_c_pw_modified_params
params_a_a      := [0.0310907, 0.01554535, 0.0168869]:
params_a_fz20   := 1.709920934161365617563962776245:
$endif

(* Equation (10) *)
g_aux := (k, rs) -> params_a_beta1[k]*sqrt(rs) + params_a_beta2[k]*rs
  + params_a_beta3[k]*rs^1.5 + params_a_beta4[k]*rs^(params_a_pp[k] + 1):
g     := (k, rs) -> -2*params_a_a[k]*(1 + params_a_alpha1[k]*rs)
  * log(1 +  1/(2*params_a_a[k]*g_aux(k, rs))):

(* Equation (8) *)
(* Attention, the function g parametrizes -alpha *)
f_pw := (rs, zeta) ->
  g(1, rs) + zeta^4*f_zeta(zeta)*(g(2, rs) - g(1, rs) + g(3, rs)/params_a_fz20)
  - f_zeta(zeta)*g(3, rs)/params_a_fz20:

f := (rs, zeta) -> f_pw(rs, zeta):
