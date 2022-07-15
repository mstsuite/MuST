(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_vxc *)

tih_par := [
  -1.0953, -0.0334789, 0.414661, -0.152399, -0.354691,
  0.0390837, -0.0748531, -0.136598, -1.41063, 0.00496577,
  0.48315, 4.02905, -0.420166, 0.0104352, -1.47409,
  -0.442455, 0.625039, 1.30351, 1.37026, -1.29598,
  -1.04305, -0.909651, -0.991782, -0.915745, 1.95026
]:

tih_zj := (j, n) -> tanh(tih_par[2*j-1] + tih_par[2*j]*n):

tih_vxc := n ->
  tih_par[17] + add(tih_par[i]*tih_zj(i-17, n), i=18..25):

f  := (rs, z) ->
  tih_vxc(n_total(rs)):
