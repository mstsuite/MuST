(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

sg4_mu1     := 0.042:
sg4_mu2     := 0.26 - sg4_mu1:
sg4_nu_MGE4 := -0.195:
sg4_k2      := -sg4_mu2^2/sg4_nu_MGE4:
sg4_k1      := 0.804 - sg4_k2:

sg4_f0 := s -> 1 + sg4_k1 + sg4_k2
   - sg4_k1*(1 - sg4_mu1*s^2/sg4_k1)/(1 - (sg4_mu1*s^2/sg4_k1)^5)
   - sg4_k2/(1 + sg4_mu2*s^2/sg4_k2):
sg4_f:= x -> sg4_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(sg4_f, rs, zeta, xs0, xs1):
