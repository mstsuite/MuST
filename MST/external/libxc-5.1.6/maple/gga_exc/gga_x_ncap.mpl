(*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_ncap_params *params;

  assert(p->params != NULL);
  params = (gga_x_ncap_params * )(p->params);
*)

ncap_f0 := s -> 1 + params_a_mu*tanh(s)*arcsinh(s)*( 1 + params_a_alpha*((1-params_a_zeta)*s*log(1+s) + params_a_zeta*s))/(1 + params_a_beta*tanh(s)*arcsinh(s)):
ncap_f  := x -> ncap_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(ncap_f, rs, zeta, xs0, xs1):
