(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_gx_params *params;

  assert(p->params != NULL);
  params = (mgga_x_gx_params * ) (p->params);
*)

$ifdef mgga_x_gx_params
params_a_c0       :=  0.827411:
params_a_c1       := -0.643560:
params_a_alphainf :=  0.852:
$endif

gx_alpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:

gx_cx0 := 4/3*(2/Pi)^(1/3):
gx_cx1 := X_FACTOR_C:

gx_gx0 := a ->
  + gx_cx0/gx_cx1
  + a*(params_a_c0 + params_a_c1*a)/(1.0 + (params_a_c0 + params_a_c1 - 1)*a) * (1 - gx_cx0/gx_cx1):

gx_gx1 := a ->
  1 + (1 - params_a_alphainf)*(1 - a)/(1 + a):

gx_f_a := a->
  + gx_gx0(a)*Heaviside(1 - a)
  + gx_gx1(a)*Heaviside(a - 1):

gx_f := (x, u, t) ->
  gx_f_a(gx_alpha(x, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(gx_f, rs, z, xs0, xs1, u0, u1, t0, t1):
