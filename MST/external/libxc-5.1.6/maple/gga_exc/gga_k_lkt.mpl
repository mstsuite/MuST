(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_k_lkt_params *params;

  assert(p->params != NULL);
  params = (gga_k_lkt_params * )(p->params);
*)

(* The m_min avoids divisions by zero *)
lkt_f0 := s -> 1/cosh(params_a_a * m_min(200, s)) + 5*s^2/3 :
lkt_f := x -> lkt_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_kinetic(lkt_f, rs, z, xs0, xs1):
