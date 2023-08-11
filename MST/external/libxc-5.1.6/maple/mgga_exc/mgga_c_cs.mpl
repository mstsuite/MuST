(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

cs_a := -0.04918:
cs_b :=  0.132:
cs_c :=  0.2533/RS_FACTOR:
cs_d :=  0.349/RS_FACTOR:

cs_thf := (z, u, t) ->
  opz_pow_n(z,8/3)*2^(-8/3)*(t - u/8):

(* This is Equation (15) of Lee1988_785 *)
(* Note that gamma = (1 - z^2) *)
f_cs := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  cs_a*(1 - z^2)/(1 + cs_d*rs) * (1 + 2*cs_b*exp(-cs_c*rs)*(
    cs_thf(z, u0, t0) + cs_thf(-z, u1, t1) - t_vw(z, xt, u0, u1)
  )):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  f_cs(rs, z, xt, xs0, xs1, u0, u1, t0, t1):
