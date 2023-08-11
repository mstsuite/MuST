(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$include "lda_c_rc04.mpl"

msigma := 1.43:
malpha := 2.30:

Bs := s -> 1/(1 + msigma*s^malpha):

f_tcs := (rs, z, xt) -> f_rc04(rs, z)*Bs(X2S*2^(1/3)*xt):

f := (rs, z, xt, xs0, xs1) ->
  f_tcs(rs, z, xt):
