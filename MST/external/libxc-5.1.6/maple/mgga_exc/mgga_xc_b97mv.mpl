(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

b97mv_par_n := 5:

b97mv_gamma_x := 0.004:
b97mv_par_x := [
    [  1.000, 0, 0],
    [  1.308, 0, 1],
    [  1.901, 0, 2],
    [  0.416, 1, 0],
    [  3.070, 1, 1]
]:

b97mv_gamma_ss := 0.2:
b97mv_par_ss := [
    [  1.000, 0, 0],
    [ -1.855, 0, 2],
    [ -5.668, 1, 0],
    [-20.497, 3, 2],
    [-20.364, 4, 2]
]:

b97mv_gamma_os := 0.006:
b97mv_par_os := [
    [  1.000, 0, 0],
    [  1.573, 0, 1],
    [ -6.298, 0, 3],
    [  2.535, 1, 0],
    [ -6.427, 3, 2]
]:

$define lda_x_params
$include "lda_x.mpl"
$include "b97mv.mpl"

b97mv_f_aux := (rs, z, xs0, xs1, ts0, ts1) ->
  + opz_pow_n( z,1)/2 * f_lda_x(rs*(2/(1 + z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_par_x,  b97mv_par_n, xs0, ts0, 0)
  + opz_pow_n(-z,1)/2 * f_lda_x(rs*(2/(1 - z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_par_x,  b97mv_par_n, xs1, ts1, 0):

f :=  (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  + b97mv_f_aux(rs, z, xs0, xs1, ts0, ts1)
  + b97mv_f(rs, z, xs0, xs1, ts0, ts1):
