(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_lypr_params *params;

  assert(p->params != NULL);
  params = (gga_c_lypr_params * )(p->params);
*)

$include "gga_c_lyp.mpl"

lypr_eta := rr -> -2/(3*sqrt(Pi))*params_a_m2*params_a_omega
  * exp(-params_a_m2^2*params_a_omega^2*rr^2):

lypr_t7 := (rr, z, xt, xs0, xs1) -> 
  -rr * (1 - z^2)/4 * (
    + 7/6*(xt^2 - lyp_aux6*(xs0^2*opz_pow_n(z,8/3) + xs1^2*opz_pow_n(-z,8/3)))
    + (1 + (1 + z)/6)*xs0^2*lyp_aux6*opz_pow_n( z, 8/3)
    + (1 + (1 - z)/6)*xs1^2*lyp_aux6*opz_pow_n(-z, 8/3)
  ):

(* This functional is very similar to gga_c_lyp. One adds the two erfc and
  the extra term proportinal to eta *)
f_lypr_rr := (rr, z, xt, xs0, xs1) -> params_a_a*(
  + erfc(params_a_m1*params_a_omega*rr)*lyp_t1(rr, z)
  + erfc(params_a_m2*params_a_omega*rr)*lyp_omega(rr)*(
    + lyp_t2(rr, z, xt) + lyp_t3(z) + lyp_t4(rr, z, xs0, xs1)
    + lyp_t5(rr, z, xs0, xs1) + lyp_t6(z, xs0, xs1)
  )
  + lyp_omega(rr)*lypr_eta(rr)*lypr_t7(rr, z, xt, xs0, xs1)
):

(* rr = rs/RS_FACTOR is equal to n_total(rs)^(-1/3) *)
f_lypr := (rs, z, xt, xs0, xs1) -> f_lypr_rr(rs/RS_FACTOR, z, xt, xs0, xs1):

f  := (rs, z, xt, xs0, xs1) -> f_lypr(rs, z, xt, xs0, xs1):

