(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_ms_params *params;

  assert(p->params != NULL);
  params = (mgga_x_ms_params * ) (p->params);
*)

ms_fa := a -> (1 - a^2)^3 / (1 + a^3 + params_a_b*a^6):
ms_f0 := (p, c) -> 1 + params_a_kappa*(1 - params_a_kappa/(params_a_kappa + MU_GE*p + c)):

ms_alpha := (t,x) -> (t - x^2/8)/K_FACTOR_C:

ms_f := (x, u, t) -> ms_f0(X2S^2*x^2, 0) + \
  ms_fa(ms_alpha(t,x))*(ms_f0(X2S^2*x^2, params_a_c) - ms_f0(X2S^2*x^2, 0)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(ms_f, rs, z, xs0, xs1, u0, u1, t0, t1):
