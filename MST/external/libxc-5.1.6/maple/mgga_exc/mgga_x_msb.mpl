(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_msb_params *params;

  assert(p->params != NULL);
  params = (mgga_x_msb_params * ) (p->params);
*)

$include "mgga_x_ms.mpl"

(* eq (5) in the paper *)
msb_beta := (t, x) -> ms_alpha(t, x)*K_FACTOR_C/(t + K_FACTOR_C):

(* eq (14) in the supplement *)
msb_fa := b -> (1 - (2*b)^2)^3 / (1 + (2*b)^3 + params_a_b*(2*b)^6):

msb_f := (x, u, t) -> ms_f0(X2S^2*x^2, 0) + \
  msb_fa(msb_beta(t,x))*(ms_f0(X2S^2*x^2, params_a_c) - ms_f0(X2S^2*x^2, 0)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(msb_f, rs, z, xs0, xs1, u0, u1, t0, t1):
