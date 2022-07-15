(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_cap_params *params;

  assert(p->params != NULL);
  params = (gga_x_cap_params * )(p->params);
*)

cap_f0 := s -> 1 - params_a_alphaoAx*s*log(1 + s)/(1 + params_a_c*log(1 + s)):
cap_f  := x -> cap_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(cap_f, rs, z, xs0, xs1):

