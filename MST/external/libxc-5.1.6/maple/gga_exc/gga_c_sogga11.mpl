(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_sogga11_params *params;

  assert(p->params != NULL);
  params = (gga_c_sogga11_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

mbeta  := 15.75592*0.004235: (* the usual value of 0.066726 *)
malpha := mbeta/(16*2^(2/3)):

sogga11_yy := (rs, z, xt) -> -malpha*mphi(z)*xt^2/(rs*f_pw(rs, z)):

sogga11_f0 := (rs, z, xt) -> 1 - 1/(1 + sogga11_yy(rs, z, xt)):
sogga11_f1 := (rs, z, xt) -> 1 - exp(-sogga11_yy(rs, z, xt)):

sogga11_t0 := (rs, z, xt) -> add(params_a_sogga11_a[i]*sogga11_f0(rs, z, xt)^(i-1), i=1..6):
sogga11_t1 := (rs, z, xt) -> add(params_a_sogga11_b[i]*sogga11_f1(rs, z, xt)^(i-1), i=1..6):

sogga11_f := (rs, z, xt, xs0, xs1) ->
  f_pw(rs, z)*(sogga11_t0(rs, z, xt) + sogga11_t1(rs, z, xt)):

f  := (rs, z, xt, xs0, xs1) -> sogga11_f(rs, z, xt, xs0, xs1):
