(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "mgga_c_scan.mpl"

(* we need to redefine two functions *)
scan_e0_g := (rs, z, t) -> (1 + 8*A(rs, z, t)*t^2)^(-1/4)/2
       + (1 + 80*A(rs, z, t)^2*t^4)^(-1/8)/2:
scan_g_infty := s -> (1 + 4*scan_chi_infty*s^2)^(-1/4)/2
       + (1 + 80*scan_chi_infty^2*s^4)^(1/8)/2:

(* and the new parameters *)
scan_b1c :=  0.030197:
scan_b2c := -0.06623:
scan_b3c :=  0.16672:

(* set parameters of f_alpha *)
params_a_c1 := 1.131:
params_a_c2 := 1.7:
params_a_d  := 1.37:
