(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

hcth_a_beta  :=  0.0042:
hcth_a_gamma :=  6:
hcth_a_c0    :=  1.09878:
hcth_a_c1    := -2.51173:
hcth_a_c2    :=  0.0156233:

hcth_a_aux := x -> 1 + hcth_a_gamma*hcth_a_beta*x*arcsinh(x):

hcth_a_f := x -> hcth_a_c0 + hcth_a_beta/X_FACTOR_C*x^2*(hcth_a_c1/hcth_a_aux(x)
     + hcth_a_c2/(hcth_a_beta*hcth_a_aux(x)^2)):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(hcth_a_f, rs, z, xs0, xs1):
