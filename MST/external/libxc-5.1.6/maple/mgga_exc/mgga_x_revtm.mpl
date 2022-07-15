(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "mgga_x_tm.mpl"

(* revtm uses the tpss definition of qtilde *)
tm_b := 0.4:
tm_qtilde := (x, t) ->
     9/20 * (tm_alpha(x, t) - 1)/sqrt(1 + tm_b*tm_alpha(x, t)*(tm_alpha(x, t) - 1))
     + 2*tm_p(x)/3:
