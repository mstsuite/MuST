(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_zvpbeint_params *params;

  assert(p->params != NULL);
  params = (gga_c_zvpbeint_params * )(p->params);
*)

params_a_gamma := (1 - log(2))/Pi^2:
params_a_BB    := 1:
$include "gga_c_pbe.mpl"

zvpbeint_nu := (rs, z, t) ->
  t*mphi(z)*(3/rs)^(1/6):

(* we write (z^2)^(omega/2) instead of z^omega in order to
   avoid the use of abs(z). Max is required not to get float
   exceptions for z->0 *)
zvpbeint_ff := (rs, z, t) ->
  exp(-params_a_alpha*zvpbeint_nu(rs, z, t)^3*m_max(z^2, 1e-20)^(params_a_omega/2)):

f  := (rs, z, xt, xs0, xs1) ->
  f_pw(rs, z) + zvpbeint_ff(rs, z, tp(rs, z, xt))*fH(rs, z, tp(rs, z, xt)):
