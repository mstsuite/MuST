(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_hl_params *params;

  assert(p->params != NULL);
  params = (lda_c_hl_params * )(p->params);
*)

$ifdef lda_c_vbh_params
params_a_hl_r := [30, 75.0]:
params_a_hl_c := [0.0252, 0.0127]:
$endif

hl_xx := (k, rs) -> rs/params_a_hl_r[k]:
hl_f0 := (k, rs) -> -params_a_hl_c[k]*
  ((1 + hl_xx(k, rs)^3)*log(1 + 1/hl_xx(k, rs)) - hl_xx(k, rs)^2 + 1/2*hl_xx(k, rs) - 1/3):

hl_f := (rs, zeta) -> hl_f0(1, rs) + f_zeta(zeta)*(hl_f0(2, rs) - hl_f0(1, rs)):
f    := (rs, zeta) -> hl_f(rs, zeta):
