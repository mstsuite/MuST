(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

a  := [0.4581652932831429, 2.217058676663745,  0.7405551735357053, 0.01968227878617998 ]:
ap := [0.119086804055547,  0.6157402568883345, 0.1574201515892867, 0.003532336663397157]:
b  := [1.0000000000000000, 4.504130959426697,  1.110667363742916,  0.02359291751427506 ]:
bp := [0.000000000000000,  0.2673612973836267, 0.2052004607777787, 0.004200005045691381]:

f  := (rs, zeta) ->
  - add((a[i] + f_zeta(zeta)*ap[i])*rs^(i-1), i=1..4) /
    add((b[i] + f_zeta(zeta)*bp[i])*rs^i,     i=1..4):
