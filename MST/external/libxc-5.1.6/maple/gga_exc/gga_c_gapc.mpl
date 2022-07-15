(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

(* Parameters from Table 6 *)
gap_par0 := [
   0.04953, 1.07924, 0.07928, (* a1, a2, a3 *)
  -2.504e-2, 7.026e-3, -1.268e-3, 1.136e-4, -3.841e-6, (* b3, b4, b5, b6, b7 *)
   0.031091,  (* params_a_a of lda_c_pw *)
   0.23878    (* pre-factor of C = 0.06483*((9*Pi)/4)^(2/3) *)
]:
gap_par1 := [
   0.0471985, 1.49676, 0.00179054,
  -3.24091e-2, 9.99978e-3, -1.93483e-3, 1.79118e-4, -6.15798e-6,
   0.015545,
   0.064535
]:

(* Equation (20): e'(rs) *)
gap_eps_1 := (rs, par) ->
  par[1]*rs^(3/2)/(1 + sqrt(rs)*(par[2] + par[3]*sqrt(rs) + par[1]*rs)):
(* Equation (21): e''(rs) *)
gap_eps_2 := (rs, par) ->
  add(par[i+1]*rs^i, i=3..7):

(* Equation (19) *)
gap_C := (rs, par) -> par[10]/rs^2:

(* Equation (17) *)
gap_c2 := (rs, z, par) ->
  + (2*f_pw(rs, z)*gap_eps_1(rs, par) - gap_C(rs, par)*gap_eps_2(rs, par))
  / (2*(gap_C(rs, par)*gap_eps_1(rs, par) - f_pw(rs, z)^2)):
(* Equation (18) *)
gap_c3 := (rs, z, par) ->
  - (2*gap_eps_1(rs, par)^2 - f_pw(rs, z)*gap_eps_2(rs, par))
  / (2*(gap_C(rs, par)*gap_eps_1(rs, par) - f_pw(rs, z)^2)):
(* Equation (16) *)
gap_c1 := (rs, z, par) ->
  - gap_C(rs, par) * gap_c3(rs, z, par):

(* after Equation (6): a = 30 is a parameter fixed by minimizing the
variance of the correlation energy error for the noble gas atoms He,
Ne, and Ar *)
gap_par_a := 30:

(* Equation (6) *)
gap_H := (rs, t, par) ->
  (gap_par_a  + par[9]*rs*log(rs)*t^2/beta_Hu_Langreth(rs))/(gap_par_a + t^2):

gap_t := (rs, z, xt) ->
  xt*n_total(rs)^(1/6)/(4*mphi(z)*(3/Pi)^(1/6)):

(* Gap function, Equation (5) *)
gap_G := (rs, z, xt, par) ->
  + mphi(z)^3*beta_Hu_Langreth(rs)*gap_t(rs, z, xt)^2
  * gap_H(rs, gap_t(rs, z, xt), par)
  / (gap_c1(rs, z, par) - gap_c2(rs, z, par)*f_pw(rs, z)):

(* Correlation energy per particle, Equation (4) *)
gap_eps := (rs, z, xt, par) ->
  + (f_pw(rs, z) + gap_c1(rs, z, par)*gap_G(rs, z, xt, par))
  / (1 + gap_c2(rs, z, par)*gap_G(rs, z, xt, par) + gap_c3(rs, z, par)*gap_G(rs, z, xt, par)^2):

(* Total energy, Equation (2) *)
f_gap := (rs, z, xt) ->
  + gap_eps(rs, 0, xt, gap_par0)
  + f_zeta(z)*(gap_eps(rs, 1, xt, gap_par1) - gap_eps(rs, 0, xt, gap_par0)):

f  := (rs, z, xt, xs0, xs1) ->
  f_gap(rs, z, xt):
