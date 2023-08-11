(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* Equations are from the Appendix of Kurth1999_889 *)

$include "gga_c_gapc.mpl"

gap_par0[10] = 0.06483*((9*Pi)/4)^(2/3): (* this is approximately equal to 0.23878 *)

(* Eq. (A4) *)
kcis_G := (rs, xt) -> xt^2*n_total(rs)^(2/3)/8:

(* Eq. (A9) *)
kcis_t := (rs, xt) -> 2^(2/3)*xt/(8*sqrt(rs)):

kcis_beta := 0.066725:

(* Eq. (A7) *)
kcis_gga0 := (rs, xt) -> f_pw(rs, 0)/(1 + kcis_beta*log(1 + kcis_t(rs, xt)^2/m_abs(f_pw(rs, 0)))):

(* Eq. (A8) *)
kcis_gga1 := (rs, xt) -> f_pw(rs, 1)/(1 + kcis_beta*log(1 + 2^(-1/3)*kcis_t(rs, xt)^2/m_abs(f_pw(rs, 1)))):

(* Eq. (A5) *)
(* The polarized parameters are the same as the unpolarized ones *)
(* except that c_1, c_2, c_3 are multiplied by 0.7. 1.5, and 2.59 respectively *)
kcis_eps_0 := (rs, xt) ->
  + (kcis_gga0(rs, xt) + gap_c1(rs, 0, gap_par0)*kcis_G(rs, xt))
  / (1 + gap_c2(rs, 0, gap_par0)*kcis_G(rs, xt) + gap_c3(rs, 0, gap_par0)*kcis_G(rs, xt)^2):

(* Eq. (A6) *)
kcis_eps_1 := (rs, xt) ->
  + (kcis_gga1(rs, xt) + 0.7*gap_c1(rs, 0, gap_par0)*kcis_G(rs, xt))
  / (1 + 1.5*gap_c2(rs, 0, gap_par0)*kcis_G(rs, xt) + 2.59*gap_c3(rs, 0, gap_par0)*kcis_G(rs, xt)^2):

(* Eq. (A2) *)
gap_f := (rs, z, xt) ->
  + kcis_eps_0(rs, xt)
  + f_zeta(z)*(kcis_eps_1(rs, xt) - kcis_eps_0(rs, xt)):

(* Eq. (A1) *)
kcis_f := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  + gap_f(rs, z, xt)
  - xs0^2/(8*ts0) * opz_pow_n( z,1)/2 * gap_f(rs*(2/(1 + z))^(1/3),  1, xs0)
  - xs1^2/(8*ts1) * opz_pow_n(-z,1)/2 * gap_f(rs*(2/(1 - z))^(1/3), -1, xs1):

f  := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  kcis_f(rs, z, xt, xs0, xs1, ts0, ts1):
