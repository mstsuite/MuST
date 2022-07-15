(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_lyp_params *params;

  assert(p->params != NULL);
  params = (gga_c_lyp_params * )(p->params);
*)

lyp_Cf := 3/10 * (3*Pi^2)^(2/3):

lyp_omega := rr -> params_a_b*exp(-params_a_c*rr)/(1 + params_a_d*rr):
lyp_delta := rr -> (params_a_c + params_a_d/(1 + params_a_d*rr))*rr:

lyp_aux6 := 1/2^(8/3):
lyp_aux4 := lyp_aux6/4:
lyp_aux5 := lyp_aux4/(9*2):

lyp_t1 := (rr, z) ->
  -(1 - z^2)/(1 + params_a_d*rr):
lyp_t2 := (rr, z, xt) ->
  -xt^2*((1 - z^2)*(47 - 7*lyp_delta(rr))/(4*18) - 2/3):
lyp_t3 := (z) ->
  -lyp_Cf/2*(1 - z^2)*(opz_pow_n(z,8/3) + opz_pow_n(-z,8/3)):
lyp_t4 := (rr, z, xs0, xs1) ->
  lyp_aux4*(1 - z^2)*(5/2 - lyp_delta(rr)/18)*(xs0^2*opz_pow_n(z,8/3) + xs1^2*opz_pow_n(-z,8/3)):
lyp_t5 := (rr, z, xs0, xs1) ->
  lyp_aux5*(1 - z^2)*(lyp_delta(rr) - 11)*(xs0^2*opz_pow_n(z,11/3) + xs1^2*opz_pow_n(-z,11/3)):
lyp_t6 := (z, xs0, xs1) ->
  -lyp_aux6*(2/3*(xs0^2*opz_pow_n(z,8/3) + xs1^2*opz_pow_n(-z,8/3))
  -opz_pow_n(z,2)*xs1^2*opz_pow_n(-z,8/3)/4 - opz_pow_n(-z,2)*xs0^2*opz_pow_n(z,8/3)/4):

f_lyp_rr := (rr, z, xt, xs0, xs1) -> params_a_a*(lyp_t1(rr, z) + lyp_omega(rr)*(
  + lyp_t2(rr, z, xt) + lyp_t3(z) + lyp_t4(rr, z, xs0, xs1)
  + lyp_t5(rr, z, xs0, xs1) + lyp_t6(z, xs0, xs1)
)):

(* rr = rs/RS_FACTOR is equal to n_total(rs)^(-1/3) *)
f_lyp := (rs, z, xt, xs0, xs1) -> f_lyp_rr(rs/RS_FACTOR, z, xt, xs0, xs1):

f  := (rs, z, xt, xs0, xs1) -> f_lyp(rs, z, xt, xs0, xs1):
