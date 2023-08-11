(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$include "vwn.mpl"

f_vwn := (rs, z) ->
  + f_aux(A_vwn[1], b_vwn[1], c_vwn[1], x0_vwn[1], rs)*(1 - f_zeta(z))
  + f_aux(A_vwn[2], b_vwn[2], c_vwn[2], x0_vwn[2], rs)*f_zeta(z):

f := (rs, z) -> f_vwn(rs, z):
