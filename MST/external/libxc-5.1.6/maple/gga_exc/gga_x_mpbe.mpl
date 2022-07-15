(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_mpbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_mpbe_params * )(p->params);
*)

mpbe_f0 := s -> s^2/(1 + params_a_a*s^2):
mpbe_f := x -> 1
  + params_a_c1*mpbe_f0(X2S*x)
  + params_a_c2*mpbe_f0(X2S*x)^2
  + params_a_c3*mpbe_f0(X2S*x)^3:

f := (rs, z, xt, xs0, xs1) -> gga_exchange(mpbe_f, rs, z, xs0, xs1):
