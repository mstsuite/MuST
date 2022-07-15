(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

tm_lambda := 0.6866:
tm_beta   := 79.873:

(* below Equation (6) *)
tm_p  := x -> (X2S*x)^2:
tm_y  := x -> (2*tm_lambda - 1)^2 * tm_p(x):

(* Equation (7) *)
tm_f0 := x -> (1 + 10*(70*tm_y(x)/27) + tm_beta*tm_y(x)^2)^(1/10):

(* after Equation (9) *)
tm_R  := (x, t) -> 1 + 595*(2*tm_lambda - 1)^2 * tm_p(x)/54 \
   - (t - 3*(tm_lambda^2 - tm_lambda + 1/2)*(t - K_FACTOR_C - x^2/72))/K_FACTOR_C:

tm_fx_DME := (x, t) -> 1/tm_f0(x)^2 + 7*tm_R(x, t)/(9*tm_f0(x)^4):

tm_alpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:

(* after Equation (11) *)
tm_qtilde := (x, t) -> 9/20*(tm_alpha(x, t) - 1) + 2*tm_p(x)/3:

(* Ratio tW/t; we have to make sure it's 1 at maximum *)
tm_tratio := (x, t) -> m_min(1.0, x^2/(8*t)):

tm_fx_SC := (x, t) -> (1 + 10*( \
       + (MU_GE + 50*tm_p(x)/729)*tm_p(x) + 146*tm_qtilde(x, t)^2/2025 \
       - 73*tm_qtilde(x,t)/405*(3/5*tm_tratio(x,t))*(1 - tm_tratio(x,t)))
       )^(1/10):

(* Equation 10 and below *)
tm_w := (x,t)-> (tm_tratio(x,t)^2 + 3*tm_tratio(x,t)^3)/(1 + tm_tratio(x,t)^3)^2:

tm_f := (x, u, t) -> tm_w(x,t)*tm_fx_DME(x, t) + (1 - tm_w(x,t))*tm_fx_SC(x, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(tm_f, rs, z, xs0, xs1, u0, u1, t0, t1):