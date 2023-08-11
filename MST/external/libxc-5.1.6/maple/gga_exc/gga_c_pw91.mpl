(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

pw91_C_c0  := 4.235e-3:
pw91_alpha := 0.09:
pw91_nu    := 16/Pi * (3*Pi^2)^(1/3):
pw91_beta  := pw91_nu*pw91_C_c0:

pw91_c1 := pw91_beta^2/(2*pw91_alpha):
pw91_c2 := 2*pw91_alpha/pw91_beta:

(* Equation (14) *)
A := (rs, z) -> pw91_c2/(exp(-2*pw91_alpha*f_pw(rs, z)/(mphi(z)^3*pw91_beta^2)) - 1):

(* Equation (13) *)
H0 := (rs, z, t) -> pw91_c1*mphi(z)^3*log(1
  + pw91_c2*(t^2 + A(rs, z)*t^4) / (1 + A(rs, z)*t^2 + A(rs, z)^2*t^4)
):

(* Pade parametrized form of C-xc found in
   M Rasolt & DJW Geldart, Phys. Rev. B 34, 1325 (1986)
*)
RS_a := [2.568, 23.266, 0.007389]:
RS_b := [1, 8.723, 0.472]:
RG_C_xc := rs -> (RS_a[1] + RS_a[2]*rs + RS_a[3]*rs^2)/(1000*(RS_b[1] + RS_b[2]*rs + RS_b[3]*rs^2)):

(* Equation (15) *)
C_xc0 := 2.568e-3:
C_x   := -0.001667:
h_a1  := -100 * 4/Pi * (4/(9*Pi))^(1/3):

H1 := (rs, z, t) -> pw91_nu * (RG_C_xc(rs) - C_xc0 - 3*C_x/7)
 * mphi(z)^3*t^2*exp(h_a1*rs*mphi(z)^4*t^2):

f_pw91 := (rs, z, xt, xs0, xs1) ->
  f_pw(rs, z) + H0(rs, z, tt(rs, z, xt)) + H1(rs, z, tt(rs, z, xt)):
f      := (rs, z, xt, xs0, xs1) -> f_pw91(rs, z, xt, xs0, xs1):
