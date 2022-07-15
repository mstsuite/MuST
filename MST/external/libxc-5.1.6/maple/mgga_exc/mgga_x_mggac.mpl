(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* replace: "mbrxc_x\(" -> "xc_mgga_x_mbrxc_get_x(" *)

$include "mgga_x_mbrxc_bg.mpl"

task_alpha := (x, t) -> (t/K_FACTOR_C) * m_max(1 - x^2/(8*t), 1e-10):

mggac_b1 := 3.712:
mggac_b2 := 2.0:
mggac_b4 := 0.1:
mggac_b3 := 2.595 + 0.5197*mggac_b4 + 0.559*mggac_b2:
mggac_b5 := -3*mggac_b3:

(* new definition of Q. The rest of the functional remains the same *)
(* We have Lambda = (32 Pi^2)^(2/3)/(6 Q) *)
mbrxc_Q := (x, t) ->
      + (32*Pi)^(2/3)/6
      * (1 + mggac_b4*task_alpha(x, t) + mggac_b5*task_alpha(x, t)^2)
      / (mggac_b1 + mggac_b2*task_alpha(x, t) + mggac_b3*task_alpha(x, t)^2):