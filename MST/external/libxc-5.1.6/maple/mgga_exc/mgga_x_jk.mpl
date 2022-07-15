(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* prefix:
  mgga_x_jk_params *params;

  assert(p->params != NULL);
  params = (mgga_x_jk_params * ) (p->params);
*)

$include "gga_x_b88.mpl"

(* equation 11 *)
y := (x,u) -> x^2 - u:
(* equation 5 *)
gBecke := x -> b88_f(x)-1:
(* equation 24 *)
jk_f := (x,u,t) -> 1 + gBecke(x)/(1 + 2*y(x,u)/x^2):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(jk_f, rs, z, xs0, xs1, u0, u1, t0, t1):