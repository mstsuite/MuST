(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

bayesian_theta0 := 1.0008:
bayesian_theta1 := 0.1926:
bayesian_theta2 := 1.8962:

bayesian_f0 := s -> s^2/(1 + s)^2:
bayesian_f  := x -> bayesian_theta0 + bayesian_f0(X2S*x)* (bayesian_theta1 + bayesian_f0(X2S*x) * bayesian_theta2):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(bayesian_f, rs, zeta, xs0, xs1):
