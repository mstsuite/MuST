(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_c_vsxc_params *params;

  assert(p->params != NULL);
  params = (mgga_c_vsxc_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

$include "gvt4.mpl"

vsxc_comp := (rs, z, spin, xs, ts) ->
  + lda_stoll_par(f_pw, rs,  z,  1)
  * gtv4(params_a_alpha_ss, params_a_dss, xs, 2*(ts - K_FACTOR_C))
  * Fermi_D(xs, ts):

(* The parallel and perpendicular components of the energy *)
vsxc_fpar  := (rs, z, xs0, xs1, ts0, ts1) ->
  + vsxc_comp(rs,  z,  1, xs0, ts0)
  + vsxc_comp(rs, -z, -1, xs1, ts1):

vsxc_fperp := (rs, z, xs0, xs1, ts0, ts1) ->
  + lda_stoll_perp(f_pw, rs,  z)
  * gtv4(params_a_alpha_ab, params_a_dab, sqrt(xs0^2 + xs1^2), 2*(ts0 + ts1 - 2*K_FACTOR_C)):

vsxc_f := (rs, z, xs0, xs1, ts0, ts1) ->
  + vsxc_fpar (rs, z, xs0, xs1, ts0, ts1)
  + vsxc_fperp(rs, z, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  vsxc_f(rs, z, xs0, xs1, ts0, ts1):
