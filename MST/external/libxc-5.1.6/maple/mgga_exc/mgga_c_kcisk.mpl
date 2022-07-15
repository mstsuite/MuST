(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "mgga_c_kcis.mpl"

kcisk_p := 0.193:
kcisk_gcnst := 20/(3*Pi*(6*Pi^2)^(1/3)):

(* Eqs. (27) *)
(* There was a "bug" in the reference implementation of Stefan Kurth. The
   f_pw was calculated for an unpolarized density instead of a polarized one.
   Here, we follow the original code. *)
kcisk_gamma0 := rs -> m_max(0, kcis_beta + 2^(1/3)*kcisk_gcnst*f_pw(2^(1/3)*rs, 0)/n_total(rs)^(1/3)):

kcisk_gamma1 := rs -> m_max(0, kcis_beta + kcisk_gcnst*f_pw(rs, 0)/n_total(rs)^(1/3)):

(* Eq. (19) and (22) *)
(* There seems to be a misspel in Eq. (22) regarding the "1 +". This follows the
   reference implementation of Stefan Kurth *)
kcis_gga0 := (rs, xt) -> f_pw(rs, 0) /
  (1 + kcisk_p*log(1 + kcisk_gamma0(rs)*kcis_t(rs, xt)^2/(kcisk_p*m_abs(f_pw(rs, 0))) )):

kcis_gga1 := (rs, xt) -> f_pw(rs, 1) /
  (1 + kcisk_p*log(1 + 2^(-1/3)*kcisk_gamma1(rs)*kcis_t(rs, xt)^2/(kcisk_p*m_abs(f_pw(rs, 1))) )):
