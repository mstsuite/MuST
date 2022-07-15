(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

c09x_mu    := 0.0617:
c09x_kappa := 1.245:
c09x_alpha := 0.0483:

c09x_f0 := s -> 1 + c09x_mu*s^2*exp(-c09x_alpha*s^2) + c09x_kappa*(1 - exp(-1/2*c09x_alpha*s^2)):
c09x_f  := x -> c09x_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(c09x_f, rs, z, xs0, xs1):
