(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* This is the definition in the paper *)
th_f0 := (x, u, t) -> -27*Pi/(10*t) * (1 + 7*x^2/(108*t)):

(* Since we write this as an enhancement functional, we need to divide
   out the LDA prefactor. The paper also defines tau without one half *)
th_f := (x, u, t) -> -th_f0(x,u,2*t) / X_FACTOR_C:

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(th_f, rs, z, xs0, xs1, u0, u1, t0, t1):
