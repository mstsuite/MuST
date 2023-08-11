(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

htbs_s1 := 0.6:
htbs_s2 := 2.6:

(* The equations to solve in order to obtain the coeficients cc are
  G(s1) = 0
  G(s2) = 1
 G'(s1) = 0
 G'(s2) = 0
G''(s1) = 0
G''(s2) = 0
*)

htbs_cc0 :=  htbs_s1^3*(htbs_s1^2 - 5*htbs_s1*htbs_s2 + 10*htbs_s2^2)/(htbs_s1 - htbs_s2)^5:
htbs_cc1 := -30*htbs_s1^2*htbs_s2^2/(htbs_s1 - htbs_s2)^5:
htbs_cc2 :=  30*htbs_s1*htbs_s2*(htbs_s1 + htbs_s2)/(htbs_s1 - htbs_s2)^5:
htbs_cc3 := -10*(htbs_s1^2 + 4*htbs_s1*htbs_s2 + htbs_s2^2)/(htbs_s1 - htbs_s2)^5:
htbs_cc4 :=  15*(htbs_s1 + htbs_s2)/(htbs_s1 - htbs_s2)^5:
htbs_cc5 := -6/(htbs_s1 - htbs_s2)^5:

$define gga_x_rpbe_params
$include "gga_x_rpbe.mpl"
$include "gga_x_wc.mpl"

htbs_g := s -> htbs_cc0 + htbs_cc1*s + htbs_cc2*s^2 + htbs_cc3*s^3 + htbs_cc4*s^4 + htbs_cc5*s^5:

htbs_f0 := s -> my_piecewise3(
   s < htbs_s1, wc_f0(s),
   my_piecewise3(s > htbs_s2, rpbe_f0(s), htbs_g(s)*rpbe_f0(s) + (1 - htbs_g(s))*wc_f0(s))
):

htbs_f  := x -> htbs_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(htbs_f, rs, z, xs0, xs1):
