(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_rtpss_params *params;

  assert(p->params != NULL);
  params = (mgga_x_rtpss_params * ) (p->params);
*)

(* These are used within the tpss_x routine *)
tpss_ff     := z -> 2:
tpss_kappa := (x, t) -> params_a_kappa:

$include "tpss_x.mpl"

(* Equation (6) *)

rtpss_f := (x, u, t) -> 1 + tpss_kappa(x, t)*(1 - exp(-tpss_fx(x, t)/tpss_kappa(x,t))):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(rtpss_f, rs, z, xs0, xs1, u0, u1, t0, t1):