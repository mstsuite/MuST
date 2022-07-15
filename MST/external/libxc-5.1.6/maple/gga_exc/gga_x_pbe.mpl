(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_pbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_pbe_params * )(p->params);
*)

(* standard PBE *)
$ifdef gga_x_pbe_params
params_a_kappa := 0.8040:
params_a_mu    := 0.2195149727645171:
$endif

(* PBE_SOL *)
$ifdef gga_x_pbe_sol_params
params_a_kappa := 0.8040:
params_a_mu    := MU_GE:
$endif

$ifdef gga_x_pbe_tca_params
params_a_kappa := 1.227:
params_a_mu    := 0.2195149727645171:
$endif

pbe_f0 := s -> 1 + params_a_kappa*(1 - params_a_kappa/(params_a_kappa + params_a_mu*s^2)):
pbe_f  := x -> pbe_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(pbe_f, rs, z, xs0, xs1):

