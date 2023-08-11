(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_k_csk_params *params;

  assert(p->params != NULL);
  params = (mgga_k_csk_params * )(p->params);
*)

csk_p := x -> X2S^2*x^2:
csk_q := u -> X2S^2*u:

(* Equation (21) *)
csk_z  := (p, q) -> 20/9*q - 40/27*p:

(* Equation (22) *)
csk_f0 := (p, q, z) ->  1 + 5*p/3 + z*csk_I(z):

(*
   The function I(z) contains exp(-1/|z|^a) which is numerically
   challenging for small |z|. Because of this, we truncate the term
   close to z=0 for negative z, since I(z) -> 1 for small |z|.
   For positive z, I(z)=1 identically, because of the step function.
*)
csk_I_negz := z -> (1 - exp(-1/abs(z)^params_a_csk_a))^(1/params_a_csk_a):
csk_I_cutoff_small := (-log(DBL_EPSILON))^(-1/params_a_csk_a):
csk_I_cutoff_large := (-log(1 - DBL_EPSILON))^(-1/params_a_csk_a):

csk_I := z -> my_piecewise5(
      z < -csk_I_cutoff_large, 0,
      z > -csk_I_cutoff_small, 1,
      csk_I_negz(m_max(m_min(z, -csk_I_cutoff_small), -csk_I_cutoff_large))
  ):

csk_f := (x, u) -> 
  csk_f0(csk_p(x), csk_q(u), csk_z(csk_p(x), csk_q(u))):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(csk_f, rs, z, xs0, xs1, u0, u1):
