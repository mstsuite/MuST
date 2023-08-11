(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

pbea_mu    := 0.00361218645365094697:
pbea_alpha := 0.52:

pbea_f  := x -> 1 + KAPPA_PBE*(1 - (1 + pbea_mu*x^2/(pbea_alpha*KAPPA_PBE))^(-pbea_alpha)):


f := (rs, z, xt, xs0, xs1) -> gga_exchange(pbea_f, rs, z, xs0, xs1):
