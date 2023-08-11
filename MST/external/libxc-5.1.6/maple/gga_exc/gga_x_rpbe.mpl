(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_rpbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_rpbe_params * )(p->params);
*)

$ifdef gga_x_rpbe_params
params_a_rpbe_kappa := KAPPA_PBE:
params_a_rpbe_mu    := MU_PBE:
$endif

rpbe_f0 := s -> 1 + params_a_rpbe_kappa * (
  1 - exp(-params_a_rpbe_mu*s^2/params_a_rpbe_kappa)
):
rpbe_f  := x -> rpbe_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(rpbe_f, rs, zeta, xs0, xs1):
