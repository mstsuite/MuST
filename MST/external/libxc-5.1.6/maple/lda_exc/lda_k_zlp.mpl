(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

f_zeta_k := z -> 1/2*(opz_pow_n(z,5/3) + opz_pow_n(-z,5/3)):

c1 := 3.2372*RS_FACTOR:
c2 := 0.00196*RS_FACTOR:

f := (rs, zeta) -> c1*f_zeta_k(zeta)/rs^2
  * (1 - c2/rs*log(1 + rs/c2)):
