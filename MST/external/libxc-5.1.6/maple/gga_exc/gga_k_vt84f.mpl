(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_k_vt84f_params *params;

  assert(p->params != NULL);
  params = (gga_k_vt84f_params * ) (p->params);
*)

(* Equation (5) *)
vt84f_f0_orig := s -> 1 - params_a_mu*s^2*exp(-params_a_alpha*s^2)/(1+params_a_mu*s^2) + (1-exp(-params_a_alpha*s^4)) * (s^(-2) - 1) + 5*s^2/3:
(* Since there's a term that looks hairy at small s, do a series expansion up to s^4 *)
vt84f_f0_series := s -> eval(convert(taylor(vt84f_f0_orig(st), st = 0, 5), polynom), st=s):
(* Glue the functions together *)
vt84f_f0 := s-> my_piecewise3(s <= sqrt(DBL_EPSILON), vt84f_f0_series(s), vt84f_f0_orig(m_max(s, sqrt(DBL_EPSILON)))):

(* Convert from x to s *)
vt84f_f := x -> vt84f_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) ->
  gga_kinetic(vt84f_f, rs, z, xs0, xs1):
