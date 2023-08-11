(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

# Note that the files have to be included in this specific order
$define gga_x_pbe_params
$include "gga_x_pbe.mpl"

$include "op.mpl"

op_qab         := 2.3789:
op_enhancement := xs -> pbe_f(xs):
