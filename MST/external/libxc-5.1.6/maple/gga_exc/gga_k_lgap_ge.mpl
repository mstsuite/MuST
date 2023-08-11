(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_k_lgap_ge_params *params;

  assert(p->params != NULL);
  params = (gga_k_lgap_ge_params * ) (p->params);
*)

(* Equation (20) *)
lgap_ge_f0 := s -> 1 + add(params_a_mu[i]*s^(i), i=1..3):
lgap_ge_f := x -> lgap_ge_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) ->
  gga_kinetic(lgap_ge_f, rs, z, xs0, xs1):
