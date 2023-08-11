(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

lg93_ad  := 1e-8:
lg93_b   := 0.024974:
lg93_a2  := (lg93_ad + 0.1234)/lg93_b:
lg93_a4  := 29.790:
lg93_a6  := 22.417:
lg93_a8  := 12.119:
lg93_a10 := 1570.1:
lg93_a12 := 55.944:

lg93_f0 := s-> 1 + lg93_a2*s^2 + lg93_a4*s^4
        + lg93_a6*s^6 + lg93_a8*s^8 + lg93_a10*s^10 + lg93_a12*s^12:

lg93_f1 := s-> lg93_f0(s)^lg93_b/(1 + lg93_ad*s^2):

lg93_f  := x->lg93_f1(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(lg93_f, rs, zeta, xs0, xs1):