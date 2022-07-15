(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

$define lda_x_params
$include "lda_x.mpl"

cc06_cnst  := (3/(4*Pi))^(2/3):

cc06_alpha := -0.0007:
cc06_beta  :=  0.0080*cc06_cnst:
cc06_gamma :=  0.026 *cc06_cnst:

cc06_f := (rs, z, us0, us1) ->
  (f_lda_x(rs, z) + f_pw(rs, z))*(1 +
    (cc06_alpha + cc06_beta*u_total(z, us0, us1))/(1 + cc06_gamma*u_total(z, us0, us1))
  ):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  cc06_f(rs, z, us0, us1):
