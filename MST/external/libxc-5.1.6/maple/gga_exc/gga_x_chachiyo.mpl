(*
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

(* x has a different definition in the Chachiyo paper  *)
chachiyo_x := x -> 2/9 * (Pi/3)**(1/3) * (2**(-1/3) * x):

(* equation 1 *)
chachiyo_f0 := x -> (3*x^2 + Pi^2*log(x+1)) / ((3*x + Pi^2)*log(x+1)):

chachiyo_f := x -> chachiyo_f0(chachiyo_x(x)):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(chachiyo_f, rs, z, xs0, xs1):
