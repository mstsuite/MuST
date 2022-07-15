(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_ev93_params *params;

  assert(p->params != NULL);
  params = (gga_x_ev93_params * )(p->params);
*)

ev93_f0 := s -> (1 + params_a_a1*s^2 + params_a_a2*s^4 + params_a_a3*s^6)/(1 + params_a_b1*s^2 + params_a_b2*s^4 + params_a_b3*s^6):
ev93_f  := x -> ev93_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(ev93_f, rs, z, xs0, xs1):
