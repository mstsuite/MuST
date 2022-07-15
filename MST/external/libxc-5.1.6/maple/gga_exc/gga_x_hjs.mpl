(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_hjs_params *params;

  assert(p->params != NULL);
  params = (gga_x_hjs_params * )(p->params);
*)

hjs_AA :=  0.757211:
hjs_BB := -0.106364:
hjs_CC := -0.118649:
hjs_DD :=  0.609650:

hjs_fH := s -> add(params_a_a[i]*s^(1+i), i=1..6)/(1 + add(params_a_b[i]*s^i, i=1..9)):

(* The m_max functions are necessary as in some cases the arguments of the
   sqrt become negative *)
hjs_zeta   := s -> m_max(s^2*hjs_fH(s), 1e-10):
hjs_eta    := s -> m_max(hjs_AA + hjs_zeta(s), 1e-10):
hjs_lambda := s -> hjs_DD + hjs_zeta(s):
hjs_chi    := (rs, z, s) -> nu(rs, z)/sqrt(hjs_lambda(s) + nu(rs, z)^2):

hjs_fF := (rs, z, s) ->
  1 - s^2/(27*hjs_CC*(1 + s^2/4)) - hjs_zeta(s)/(2*hjs_CC):

hjs_fG := (rs, z, s) ->
  - 2/5  * hjs_CC*hjs_fF(rs, z, s)*hjs_lambda(s)
  - 4/15 * hjs_BB*hjs_lambda(s)^2
  - 6/5  * hjs_AA*hjs_lambda(s)^3
  - hjs_lambda(s)^(7/2)*(4/5*sqrt(Pi) + 12/5*(sqrt(hjs_zeta(s)) - sqrt(hjs_eta(s)))):

hjs_f1 := (rs, z, s) ->
   + hjs_AA
   - 4/9 * hjs_BB*(1 - hjs_chi(rs, z, s))/hjs_lambda(s)
   - 2/9 * hjs_CC*hjs_fF(rs, z, s)*(2 - 3*hjs_chi(rs, z, s) + hjs_chi(rs, z, s)^3)/hjs_lambda(s)^2
   - 1/9 * hjs_fG(rs, z, s)*(8 - 15*hjs_chi(rs, z, s) + 10*hjs_chi(rs, z, s)^3 - 3*hjs_chi(rs, z, s)^5)/hjs_lambda(s)^3
   + 2*nu(rs, z)*(sqrt(hjs_zeta(s) + nu(rs, z)^2) -  sqrt(hjs_eta(s) + nu(rs, z)^2))
   + 2*hjs_zeta(s)*log((nu(rs, z) + sqrt(hjs_zeta(s) + nu(rs, z)^2))/(nu(rs, z) + sqrt(hjs_lambda(s) + nu(rs, z)^2)))
   - 2*hjs_eta(s)*log((nu(rs, z) + sqrt(hjs_eta(s) + nu(rs, z)^2))/(nu(rs, z) + sqrt(hjs_lambda(s) + nu(rs, z)^2))):

hjs_fx := (rs, z, x) -> hjs_f1(rs, z, X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange_nsp(hjs_fx, rs, z, xs0, xs1):
