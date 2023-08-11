(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* replace: "my_dilog\(" -> "xc_dilogarithm(" *)

`diff/my_dilog` := proc(g, x) -log(1 - g)/g * diff(g, x) end proc:

gg99_a := 3^(1/4)/(2*sqrt(2)*Pi^(3/2)):
gg99_b := 4*sqrt(3)*Pi^3:

(* Equation 22 in the paper, i.e.
   the solution of x = 2*Pi*sinh(r)/(3*cosh(r))^(1/3) *)

gg99_r_branch1 := x -> arcsinh( (gg99_a * x * sqrt(x^2 + (gg99_b + sqrt(gg99_b^2 - x^6))^(2/3))) / (gg99_b + sqrt(gg99_b^2 - x^6))^(1/6) ):

(* The second branch is from Andrew Gilbert via email *)

gg99_r_branch2 := x-> arcsinh(sqrt(x^3*(3/gg99_b)*cos(arctan(sqrt(1/(gg99_b^2)*x^6-1))/3))):

(* Glue the pieces together. The min and max are required
   to avoid float exceptions *)
gg99_r := x -> my_piecewise3(x < gg99_b^(1/3),
  gg99_r_branch1(m_min(x, gg99_b^(1/3) - 1e-10)),
  gg99_r_branch2(m_max(x, gg99_b^(1/3) + 1e-10))
  ):

(* Equation 21 *)
gg99_f0 := r -> (Pi^2 - 12*r*log(1 + exp(-2*r)) +
12*my_dilog(-exp(-2*r))) / (2*3^(1/3)*Pi*r*sech(r)^(2/3)) / X_FACTOR_C:

(* Assemble the function *)
gg99_f := x -> gg99_f0(gg99_r(x)):

f := (rs, zeta, xt, xs0, xs1) ->
  gga_exchange(gg99_f, rs, zeta, xs0, xs1):
