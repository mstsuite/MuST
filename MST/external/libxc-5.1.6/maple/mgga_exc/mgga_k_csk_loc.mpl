(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_k_csk_loc_params *params;

  assert(p->params != NULL);
  params = (mgga_k_csk_loc_params * )(p->params);
*)

$include "mgga_k_csk.mpl"

(* Equation (21) *)
csk_z  := (p, q) -> 1 + params_a_csk_cp*p + params_a_csk_cq*q - (1 + 5*p/3):