(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_pz_params *params;

  assert(p->params != NULL);
  params = (lda_c_pz_params * )(p->params);
*)

$ifdef lda_c_pz_params
params_a_gamma := [-0.1423, -0.0843]:
params_a_beta1 := [ 1.0529,  1.3981]:
params_a_beta2 := [ 0.3334,  0.2611]:
params_a_a     := [ 0.0311,  0.01555]:
params_a_b     := [-0.048,  -0.0269]:
params_a_c     := [ 0.0020,  0.0007]:
params_a_d     := [-0.0116, -0.0048]:
$endif

(* Equation C3 *)
ec_low  := (i, rs) -> params_a_gamma[i] / \
        (1 + params_a_beta1[i]*sqrt(rs) + params_a_beta2[i]*rs):

(* Equation [1].C5 *)
ec_high := (i, rs) -> params_a_a[i]*log(rs) + params_a_b[i] \
        + params_a_c[i]*rs*log(rs) + params_a_d[i]*rs:

ec := (i, x) -> my_piecewise3(x >= 1, ec_low(i, x), ec_high(i, x)):

f_pz := (rs, zeta) -> \
 ec(1, rs) + (ec(2, rs) - ec(1, rs))*f_zeta(zeta):

f := (rs, zeta) -> f_pz(rs, zeta):
