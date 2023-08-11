(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_mn12_params *params;

  assert(p->params != NULL);
  params = (mgga_x_mn12_params * ) (p->params);
*)

$define lda_x_params
$include "lda_x.mpl"

mn12_omega_x := 2.5:
mn12_gamma_x := 0.004:

mn12_vx := (rs, z) -> 1/(1 + rs/(mn12_omega_x*RS_FACTOR)*(2/(1 + z))^(1/3)):
mn12_ux := x -> mn12_gamma_x*x^2/(1 + mn12_gamma_x*x^2):
mn12_wx := t -> (K_FACTOR_C - t)/(K_FACTOR_C + t):

mn12_pol1  := t-> params_a_c[ 1] + params_a_c[ 2]*mn12_wx(t) + params_a_c[ 3]*mn12_wx(t)^2 + params_a_c[ 4]*mn12_wx(t)^3
  + params_a_c[ 5]*mn12_wx(t)^4 + params_a_c[ 6]*mn12_wx(t)^5:
mn12_pol2  := t-> params_a_c[ 7] + params_a_c[ 8]*mn12_wx(t) + params_a_c[ 9]*mn12_wx(t)^2 + params_a_c[10]*mn12_wx(t)^3
   + params_a_c[11]*mn12_wx(t)^4:
mn12_pol3  := t-> params_a_c[12] + params_a_c[13]*mn12_wx(t) + params_a_c[14]*mn12_wx(t)^2 + params_a_c[15]*mn12_wx(t)^3:
mn12_pol4  := t-> params_a_c[16] + params_a_c[17]*mn12_wx(t) + params_a_c[18]*mn12_wx(t)^2:
mn12_pol5  := t-> params_a_c[19] + params_a_c[20]*mn12_wx(t) + params_a_c[21]*mn12_wx(t)^2 + params_a_c[22]*mn12_wx(t)^3
   + params_a_c[23]*mn12_wx(t)^4:
mn12_pol6  := t-> params_a_c[24] + params_a_c[25]*mn12_wx(t) + params_a_c[26]*mn12_wx(t)^2 + params_a_c[27]*mn12_wx(t)^3:
mn12_pol7  := t-> params_a_c[28] + params_a_c[29]*mn12_wx(t) + params_a_c[30]*mn12_wx(t)^2:
mn12_pol8  := t-> params_a_c[31] + params_a_c[32]*mn12_wx(t) + params_a_c[33]*mn12_wx(t)^2 + params_a_c[34]*mn12_wx(t)^3:
mn12_pol9  := t-> params_a_c[35] + params_a_c[36]*mn12_wx(t) + params_a_c[37]*mn12_wx(t)^2:
mn12_pol10 := t-> params_a_c[38] + params_a_c[39]*mn12_wx(t) + params_a_c[40]*mn12_wx(t)^2:

mn12_f := (rs, z, x, u, t) ->
  + mn12_pol1(t)
  + mn12_pol2(t)*mn12_ux(x)
  + mn12_pol3(t)*mn12_ux(x)^2
  + mn12_pol4(t)*mn12_ux(x)^3
  + mn12_pol5(t)*mn12_vx(rs, z)
  + mn12_pol6(t)*mn12_ux(x)*mn12_vx(rs, z)
  + mn12_pol7(t)*mn12_ux(x)^2*mn12_vx(rs, z)
  + mn12_pol8(t)*mn12_vx(rs, z)^2
  + mn12_pol9(t)*mn12_ux(x)*mn12_vx(rs, z)^2
  + mn12_pol10(t)*mn12_vx(rs, z)^3:

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange_nsp(mn12_f, rs, z, xs0, xs1, u0, u1, t0, t1):
