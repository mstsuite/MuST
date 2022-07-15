(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_k_rational_p_params *params;

  assert(p->params != NULL);
  params = (gga_k_rational_p_params * )(p->params);
*)

rational_p_f0 := s -> (1 + params_a_C2/params_a_p * s^2)^(-params_a_p):
rational_p_f := x -> rational_p_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_kinetic(rational_p_f, rs, z, xs0, xs1):
