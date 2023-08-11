(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

op_a1 := 1.5214:
op_a2 := 0.5764:
op_b1 := 1.1284:
op_b2 := 0.3183:

(* This wrapper is to avoid overflows in the OP functionals. The
   energy is not affected, since the value is only changed for
   densities that are screened away.  *)
op_b88_zab := (f_x, rs, z, xs0, xs1) ->
  my_piecewise3(
    b88_zab(1, op_enhancement, rs, z, xs0, xs1) = 0,
    DBL_EPSILON,
    b88_zab(1, op_enhancement, rs, z, xs0, xs1)
  ):

op_beta := (rs, z, xs0, xs1) ->
  op_qab/op_b88_zab(op_enhancement, rs, z, xs0, xs1):

op_f_s := (rs, z, xt, xs0, xs1) ->
  - (1 - z^2)*n_total(rs)/4.0
  * (op_a1*op_beta(rs, z, xs0, xs1) + op_a2)
  / (op_beta(rs, z, xs0, xs1)^4 + op_b1*op_beta(rs, z, xs0, xs1)^3 + op_b2*op_beta(rs, z, xs0, xs1)^2)
:

op_f := (rs, z, xt, xs0, xs1) ->
   my_piecewise3(1 - abs(z) <= p_a_zeta_threshold or (screen_dens(rs,z) and screen_dens(rs,-z)), 0,
                 op_f_s(rs, z_thr(z), xt, xs0, xs1)):

f := (rs, z, xt, xs0, xs1) ->
  op_f(rs, z, xt, xs0, xs1):
