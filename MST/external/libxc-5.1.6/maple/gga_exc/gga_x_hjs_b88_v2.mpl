(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_hjs_params *params;

  assert(p->params != NULL);
  params = (gga_x_hjs_params * )(p->params);
*)

$include "gga_x_hjs.mpl"

hjs2_xi := 1/(exp(20) - 1):
hjs2_fs := s -> -log((exp(-s) + hjs2_xi)/(1 + hjs2_xi)):

hjs_fx := (rs, z, x) -> hjs_f1(rs, z, hjs2_fs(X2S*x)):
