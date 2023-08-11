(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_vxc *)

(* prefix:
  mgga_x_tb09_params *params;

  assert(p->params != NULL);
  params = (mgga_x_tb09_params * ) (p->params);
*)

(* replace: "br89_x\(" -> "xc_mgga_x_br89_get_x(" *)

$include "mgga_x_br89.mpl"

params_a_gamma := 0.8:
tb09_c_HEG := (3*params_a_c - 2)*sqrt(5/12)/(Pi):

# we add here a threshold of 1e-10 for either tau or the Fermi curvature
tb09_f := (rs, z, x, u, t) -> (params_a_c*X_FACTOR_C*br89_v(br89_x(br89_cQ(br89_Q(x, u, t))))
  + tb09_c_HEG*sqrt(2*m_max(t - params_a_alpha*x^2/8, 1e-10)))*n_spin(rs, z)^(1/3):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  tb09_f(rs, z, xs0, u0, t0):
