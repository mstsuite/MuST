(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

gvt4_gamm := (alpha, x, z) -> 1 + alpha*(x^2 + z):

gtv4 := (alpha, dd, x, z) ->
  dd[1]/gvt4_gamm(alpha, x, z) +
  (dd[2]*x^2 + dd[3]*z)/gvt4_gamm(alpha, x, z)^2 +
  (dd[4]*x^4 + dd[5]*x^2*z + dd[6]*z^2)/gvt4_gamm(alpha, x, z)^3:
