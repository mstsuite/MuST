(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_k_gds08_params *params;

  assert(p->params != NULL);
  params = (lda_k_gds08_params * )(p->params);
*)

gds08_fs := (rs, z) -> (1 + z)/2 *(
  + params_a_A
  + params_a_B*log(2*n_spin(rs, z))
  + params_a_C*log(2*n_spin(rs, z))^2
):

# Eq. (12)
gds08_f := (rs, z) ->
  + my_piecewise3(screen_dens(rs,  z), 0, gds08_fs(rs, z_thr( z)))
  + my_piecewise3(screen_dens(rs, -z), 0, gds08_fs(rs, z_thr(-z))):

f := (rs, z) -> gds08_f(rs, z):