(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_pbeint_params *params;

  assert(p->params != NULL);
  params = (gga_x_pbeint_params * )(p->params);
*)

pbeint_mu := s -> params_a_muGE + (params_a_muPBE - params_a_muGE)* \
   params_a_alpha*s^2/(1 + params_a_alpha * s^2):

(* this is the gga_x_pbe expression *)
pbeint_f0 := s -> 1 + params_a_kappa * (1 - params_a_kappa/(params_a_kappa + pbeint_mu(s)*s^2)):
pbeint_f  := x -> pbeint_f0(X2S * x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(pbeint_f, rs, z, xs0, xs1):
