(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

cs1_gamma :=  0.006: (* as in B88 *)
cs1_d     :=  0.349: (* as in CS  *)
cs1_C1    := -0.018897:
cs1_C2    :=  0.155240:
cs1_C3    := -0.159068:
cs1_C4    :=  0.007953:

(* Equation (24) corrected in Equation (8) in Proynov2006_436 *)
cs1_ess := (rs, z, xs) ->
  + opz_pow_n(z,1)/2 * n_spin(rs, z)^(1/3)/(n_spin(rs, z)^(1/3) + cs1_d)
  * (cs1_C1 + cs1_C2*cs1_gamma^2*xs^4/(1 + cs1_gamma*xs^2)^2):

(* Equation (25) corrected in Equation (6) in Proynov2006_436 *)
cs1_eab := (rs, z, xt) ->
  + (1 - z^2)/4 * 1/(1 + cs1_d*n_total(rs)^(-1/3))
  * (cs1_C3 + cs1_C4*cs1_gamma^2*xt^4/(1 + cs1_gamma*xt^2)^2):

f_cs1 := (rs, z, xt, xs0, xs1) ->
  + cs1_eab(rs,  z, xt)
  + cs1_ess(rs,  z, xs0)
  + cs1_ess(rs, -z, xs1):

f := (rs, z, xt, xs0, xs1) ->
  f_cs1(rs, z, xt, xs0, xs1):
