(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* replace: "BesselK\(0," -> "xc_bessel_K0(" *)
(* replace: "BesselK\(1," -> "xc_bessel_K1(" *)
(* replace: "int1\(" -> "xc_integrate(func1, NULL, 0.0, " *)
(* replace: "int2\(" -> "xc_integrate(func2, NULL, 0.0, " *)
(* prefix:
  lda_x_1d_exponential_params *params;

  assert(p->params != NULL);
  params = (lda_x_1d_exponential_params * )(p->params);
*)

$include "lda_x_1d_exponential.mpl"

x1d_inter := x -> 2.0*BesselK(0, x):
