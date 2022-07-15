(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

pearson_f0 := s -> 1 + 5/27*s^2/(1 + s^6):
pearson_f  := x -> pearson_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(pearson_f, rs, zeta, xs0, xs1):
