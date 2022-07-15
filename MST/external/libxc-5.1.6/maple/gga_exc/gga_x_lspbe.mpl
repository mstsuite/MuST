(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_lspbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_lspbe_params * )(p->params);
*)

lspbe_f0 := s -> 1 + params_a_kappa*(1 - params_a_kappa/(params_a_kappa + params_a_mu*s^2))
            - (params_a_kappa+1)*(1-exp(-params_a_alpha*s^2)):
lspbe_f  := x -> lspbe_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(lspbe_f, rs, z, xs0, xs1):

