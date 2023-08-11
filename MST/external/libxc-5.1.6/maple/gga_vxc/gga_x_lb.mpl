(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_vxc *)
(* prefix:
  gga_x_lb_params *params;

  assert(p->params != NULL);
  params = (gga_x_lb_params * )(p->params);
*)

lb_f0 := (rs, z, x) -> -my_piecewise3(x < 300,
         params_a_beta*x^2/(1 + 3*params_a_beta*x*arcsinh(params_a_gamma*x)),
         x/(3*log(2*params_a_gamma*x))):

lb_f := (rs, z, x) -> (params_a_alpha*(4/3)*LDA_X_FACTOR + lb_f0(rs, z, x))*n_spin(rs, z)^(1/3):

f := (rs, z, xt, xs0, xs1) -> lb_f(rs, z, xs0):
