(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* nop scaling *)
s_scaling_0 := s -> s:

(* original scaling of Heyd *)
strans :=  8.3:
smax   :=  8.572844:
sconst := 18.79622316:
s_scaling_1 := s -> my_piecewise3(
  s < strans, s,
  smax - sconst/s^2
):

(* first version of the scaling by TM Henderson, apparently used by Gaussian *)
s_scaling_2 := s -> my_piecewise3(
  s < 1,  s,
  my_piecewise3(s > 15, smax, m_max(m_min(s, 15), 1) - log(1 + exp(m_max(m_min(s, 15), 1) - smax)))
):

(* second version of the scaling by TM Henderson *)
s_scaling_3 := s -> s - (1 - exp(-s))*log(1 + exp(s - smax)):

(* appendix of JCP 128, 194105 (2008) *)
s_p :=  [0.615482, 1.136921, -0.449154, 0.0175739*8.572844]:
s_q :=  [1.229195, -0.0269253, 0.313417, -0.0508314, 0.0175739]:

s_scaling_4 := s->
  s*(1 + s^3*add(s_p[i]*s^i, i=1..4))/(1 + s^3*add(s_q[i]*s^i, i=1..5)):
