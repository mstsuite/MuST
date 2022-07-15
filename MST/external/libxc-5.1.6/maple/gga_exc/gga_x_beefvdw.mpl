(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

beefvdw_coeffs := [
   1.516501714e0,   4.413532099e-1, -9.182135241e-2, -2.352754331e-2,
   3.418828455e-2,  2.411870076e-3, -1.416381352e-2,  6.975895581e-4,
   9.859205137e-3, -6.737855051e-3, -1.573330824e-3,  5.036146253e-3,
  -2.569472453e-3, -9.874953976e-4,  2.033722895e-3, -8.018718848e-4,
  -6.688078723e-4,  1.030936331e-3, -3.673838660e-4, -4.213635394e-4,
   5.761607992e-4, -8.346503735e-5, -4.458447585e-4,  4.601290092e-4,
  -5.231775398e-6, -4.239570471e-4,  3.750190679e-4,  2.114938125e-5,
  -1.904911565e-4,  7.384362421e-5
]:

beefvdw_k  := 4:
beefvdw_xi := p -> 2*p/(beefvdw_k + p) - 1:

with(orthopoly):
beefvdw_f := x -> add(beefvdw_coeffs[i]*P(i-1, beefvdw_xi(X2S^2*x^2)), i=1..30):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(beefvdw_f, rs, zeta, xs0, xs1):
