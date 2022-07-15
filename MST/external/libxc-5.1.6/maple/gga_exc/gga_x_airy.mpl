(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

airy_a5  := 133.983631:
airy_a6  :=   3.217063:
airy_a7  := 136.707378:
airy_a8  :=   3.223476:
airy_a9  :=   2.675484:
airy_a10 :=   3.473804:

airy_f1 := s -> (1 - airy_a5*s^airy_a6 + airy_a7*s^airy_a8)/(1 + airy_a9*s^airy_a10):

$include "gga_x_lag.mpl"

airy_f := x -> lag_f0(X2S*x) + airy_f1(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(airy_f, rs, zeta, xs0, xs1):
