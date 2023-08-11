(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_ityh_optx_params *params;

  assert(p->params != NULL);
  params = (gga_x_ityh_optx_params * )(p->params);
*)

$include "gga_x_optx.mpl"
$include "gga_x_ityh.mpl"

ityh_enhancement := xs -> optx_f(xs):
