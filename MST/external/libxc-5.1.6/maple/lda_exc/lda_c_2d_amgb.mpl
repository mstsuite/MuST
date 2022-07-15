(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$define xc_dimensions_2d

amgb_aa := [ -0.1925,     0.117331,    0.0234188 ]:
amgb_bb := [  0.0863136, -3.394e-2,   -0.037093  ]:
amgb_cc := [  0.0572384, -7.66765e-3,  0.0163618 ]:
amgb_ee := [  1.0022,     0.4133,      1.424301  ]:
amgb_ff := [ -0.02069,    0,         0       ]:
amgb_gg := [  0.33997,    6.68467e-2,  0       ]:
amgb_hh := [  1.747e-2,   7.799e-4,    1.163099  ]:

amgb_ax = -4/(3*Pi*sqrt(2)):
amgb_beta := 1.3386:

amgb_dd := [seq(-amgb_aa[i]*amgb_hh[i], i=1..3)]:

amgb_alpha := (i, rs) -> amgb_aa[i]
  + (amgb_bb[i]*rs + amgb_cc[i]*rs^2 + amgb_dd[i]*rs^3)
  * log(1 + 1/(amgb_ee[i]*rs + amgb_ff[i]*rs^1.5 + amgb_gg[i]*rs^2 + amgb_hh[i]*rs^3)):

amgb_ex6 := (rs, z) -> -4*sqrt(2)/(3*Pi*rs)
  * (f_zeta_2d(z) - 1 - 3/8*z^2 - 3/128*z^4):

f_amgb := (rs, z) ->
  amgb_alpha(1, rs) + amgb_alpha(2, rs)*z^2 + amgb_alpha(3, rs)*z^4 + (exp(-amgb_beta*rs) - 1)*amgb_ex6(rs, z):

f := (rs, z) -> f_amgb(rs, z):
