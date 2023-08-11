(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_k_pg_params *params;

  assert(p->params != NULL);
  params = (gga_k_pg_params * )(p->params);
*)

pg_f0 := s -> 5/3*s^2 + exp(-params_a_pg_mu * s^2):
pg_f := x -> pg_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_kinetic(pg_f, rs, z, xs0, xs1):
