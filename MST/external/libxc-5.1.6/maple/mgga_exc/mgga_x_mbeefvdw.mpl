(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

mbeef_n := 5:
mbeef_coefs := [
    [ 1.17114923e+00,  1.15594371e-01, -5.32167416e-02, -2.01131648e-02,  1.41417107e-03],
    [-6.76157938e-02,  4.53837246e-02, -2.22650139e-02,  1.92374554e-02,  9.19317034e-07],
    [ 1.48659502e-02,  3.18024096e-02, -5.21818079e-03,  1.33707403e-07, -5.00749348e-07],
    [ 1.40794142e-03, -6.08338264e-03, -6.57949254e-07, -5.49909413e-08,  5.74317889e-08],
    [ 1.41530486e-04, -1.00478906e-07,  2.01895739e-07,  3.97324768e-09, -3.40722258e-09]
]:

$include "mbeef.mpl"

mbeefvdw_f := (x, u, t) -> mbeef_expansion(x, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(mbeefvdw_f, rs, z, xs0, xs1, u0, u1, t0, t1):
