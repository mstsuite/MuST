(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
$include "mgga_c_scan.mpl"

(* interpolation function, defined in main text *)
rregtm_gamma1 := 0.2:
rregtm_g := a -> (1 + rregtm_gamma1)*a / (rregtm_gamma1 + a):
rregtm_f2g := g -> 3*g^3 / (1 + g^3 + g^6):
rregtm_f2 := a -> rregtm_f2g(rregtm_g(a)):

(* equation 3 *)
rregtm_f := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  f_pbe(rs, z, xt, xs0, xs1) + rregtm_f2(scan_alpha(z, xt, ts0, ts1))*(
    + scan_e0(rs, z, X2S*2^(1/3)*xt)
    - f_pbe(rs, z, xt, xs0, xs1)
  ):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  scan_f(rs, z, xt, xs0, xs1, ts0, ts1):
