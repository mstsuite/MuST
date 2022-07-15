(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

(* constants from text in section 2 *)
pbetrans_kappa_pbe := 0.814: (* This is probably a misprint in the manuscript as KAPPA_PBE is 0.8040 *)
pbetrans_kappa_revpbe := 1.227:
pbetrans_mu := 0.219:

(* parameters from section 4 *)
pbetrans_alpha := 2*(3*Pi^2)^(1/3):
pbetrans_beta := 3:

(* eq 3 *)
pbetrans_fermi := s -> 1/(1+exp(-pbetrans_alpha*(s-pbetrans_beta))):
(* eq 5 *)
pbetrans_kappa := s -> (1-pbetrans_fermi(s))*pbetrans_kappa_revpbe + pbetrans_fermi(s)*pbetrans_kappa_pbe:
(* eq 4 *)
pbetrans_f0 := s -> 1 + pbetrans_kappa(s)*(1 - pbetrans_kappa(s)/(pbetrans_kappa(s) + pbetrans_mu*s^2)):

pbetrans_f  := x -> pbetrans_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(pbetrans_f, rs, z, xs0, xs1):


