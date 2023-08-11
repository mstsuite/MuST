(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

C0 := (1 - log(2))/(2*Pi^2): (* Equation (9) *)
C1 := 4*C0/3:                (* Equations (13), (28), (33) *)
C2 := RS_FACTOR:             (* Equation (8) *)
C3 := C2/3:

(* several cutoffs *)
big     := 1.0e4:
cutoff  := 1.0e7:
ei_xmax := 7.0183341467e+02:

(* Equation (39) *)
kssp0_k0 := 1.291551074:
kssp0_k1 := 0.349064173:
kssp0_r1 := 0.08327588:

kssp0 := rs ->
  kssp0_k0 - kssp0_k1*(1 - exp(-kssp0_r1*rs^(4/5))):

(* Equation (45) *)
fssp_A1 := 1.622118767:
fssp_A2 := 0.489958076:
fssp_A3 := 1.379021941:

fssp := (rs, gr) ->
  (1 + fssp_A1*gr + fssp_A2^2*gr^2)*exp(-fssp_A2^2*gr^2)/sqrt(1 + fssp_A3*gr/rs):

(* Equation (34) *)
fa_a1 := 0.939016:
fa_a2 := 1.733170:

f_factor := rs -> exp(-rs^2/(fa_a1*sqrt(rs) + fa_a2*rs)^2):

(* Equation (40) *)
kss0_k0 :=  1.200801774:
kss0_k1 :=  0.859614445:
kss0_k2 := -0.812904345:
kss0_r1 :=  1.089338848:
kss0_r2 :=  0.655638823:

kss0 := (rs, gr) ->
  + kss0_k0
  + kss0_k1*(1 - exp(-kss0_r1*sqrt(rs)))
  + kss0_k2*(1 - exp(-kss0_r2*rs^(2/5))):

fss_A4 := 4.946281353:
fss_A5 := 3.600612059:

fss := (rs, gr) ->
  (1 + fss_A4^2*gr^2)*exp(-fss_A4^2*gr^2)/sqrt(1 + fss_A5*gr/rs):

(* Equation (15) *)
eq15 := mu -> (3 + 2*(sqrt(mu) + mu))/(3 + 6*(sqrt(mu) + mu)):

f_eab := mu ->
  C0*(-xc_E1_scaled(mu)*(1 + 2*mu*eq15(mu)) + 2*eq15(mu)):

(* Equation (13) *)
(*
   This is numerically suboptimal - the max function cuts off f_eab to
   around 10^-6. This is too large. The other possibility that was used
   before in the code was to put term to zero if mu > ei_xmax
*)
mu_ba := (rsa, ga2) -> C1*rsa/m_max(kssp0(rsa)^2*fssp(rsa, ga2)^2, 1e-60):
term1 := (rsa, z, ga2) -> my_piecewise3(mu_ba(rsa, ga2) >= cutoff, 0,
  f_eab(mu_ba(rsa, ga2))*(1 - z)/2):


mu_aa := (rsa, ga2) -> C1*rsa/m_max(kss0(rsa)^2*fss(rsa, ga2)^2, 1e-60):
term2 := (rsa, z, ga2) -> my_piecewise3(mu_aa(rsa, ga2) >= cutoff, 0,
  f_eab(mu_aa(rsa, ga2))*f_factor(rsa)*(1 + z)/2):

f_ft97 := (rs, z, xs) ->
  + term1(rs*(2/(1 + z))^(1/3), z, C3^2*xs^2)
  + term2(rs*(2/(1 + z))^(1/3), z, C3^2*xs^2):

f  := (rs, z, xt, xs0, xs1) -> f_ft97(rs, z, xs0) + f_ft97(rs, -z, xs1):
