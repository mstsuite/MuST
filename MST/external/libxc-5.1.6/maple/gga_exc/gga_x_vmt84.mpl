(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_vmt84_params *params;

  assert(p->params != NULL);
  params = (gga_x_vmt84_params * )(p->params);
*)

$include "gga_x_vmt.mpl"

vmt84_f0 := s -> (1 - exp(-params_a_alpha*s^4))/s^2 - 1 + exp(-params_a_alpha*s^4):
vmt84_f  := x -> vmt_f(x) + vmt84_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(vmt84_f, rs, zeta, xs0, xs1):
