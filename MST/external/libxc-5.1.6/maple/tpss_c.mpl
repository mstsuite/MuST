(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* Equation 28 squared with Equation 25 built in *)
tpss_xi2 := (z, xt, xs0, xs1) ->
  (1 - z^2)*(t_total(z, xs0^2, xs1^2) - xt^2)/(2*(3*Pi^2)^(1/3))^2:

(* Equation 33 *)
tpss_C00 := (cc, z) ->
  + add(cc[i]*z^(2*(i-1)), i=1..4):

(* Equation 34 *)
(* The series expansion of C0 for z -> +- 1 goes as C + D*(1+-z)^(2/3) + O(1+-z)
   whose first derivative seems to diverge for ferromagnetic densities, leading
   to severe numerical instabilities. I can not correct for the bad design of
   the functional, and this is the possible solution that I found *)
tpss_C0_den := (z, xt, xs0, xs1) ->
  1 + tpss_xi2(z, xt, xs0, xs1)*((1 + z)^(-4/3) + (1 - z)^(-4/3))/2:
tpss_C0 := (cc, z, xt, xs0, xs1) -> my_piecewise3(1 - m_abs(z) <= 1e-12,
  add(cc[i], i=1..4),
  tpss_C00(cc, z) / tpss_C0_den(z_thr(z), xt, xs0, xs1)^4):

(* Equation 11, with tau_W from Equation 12 *)
tpss_aux := (z, xt, ts0, ts1) ->
  m_min(xt^2/(8*t_total(z, ts0, ts1)), 1):

(* n_sigma/n \epsilon^sigma in Equation 25 *)
tpss_par_s0 := (f_gga, rs, z, xt, xs0, xs1) ->
  m_max(f_gga(rs*(2/(1 + z))^(1/3),  1, xs0, xs0, 0), f_gga(rs, z, xt, xs0, xs1))*(1 + z)/2:
tpss_par_s1 := (f_gga, rs, z, xt, xs0, xs1) ->
  m_max(f_gga(rs*(2/(1 - z))^(1/3), -1, xs1, 0, xs1), f_gga(rs, z, xt, xs0, xs1))*(1 - z)/2:

(* Second line of Equation 25 *)
(* The screening of the density is important in order to stabilize this functional
   for ferromagnetic densities *)
tpss_par  := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  - (1 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1))*tpss_aux(z, xt, ts0, ts1)^2*(
    + my_piecewise3(screen_dens_zeta(rs,  z),
        f_gga(rs, z_thr(z), xt, xs0, xs1)*(1 + z)/2,
        tpss_par_s0(f_gga, rs, z_thr(z), xt, xs0, xs1)
      )
    + my_piecewise3(screen_dens_zeta(rs, -z),
        f_gga(rs, z_thr(z), xt, xs0, xs1)*(1 - z)/2,
        tpss_par_s1(f_gga, rs, z_thr(z), xt, xs0, xs1)
      )
  ):

(* First line of Equation 25 *)
tpss_perp := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  (1 + tpss_C0(params_a_C0_c, z, xt, xs0, xs1)*tpss_aux(z, xt, ts0, ts1)^2)
  * f_gga(rs, z, xt, xs0, xs1):

(* Equation in full 25 *)
tpss_f0 := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + tpss_par (f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  + tpss_perp(f_gga, rs, z, xt, xs0, xs1, ts0, ts1):

(* Equation 24 *)
tpss_f := (f_gga, rs, z, xt, xs0, xs1, ts0, ts1) ->
  + tpss_f0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)
  * (1 + params_a_d*tpss_f0(f_gga, rs, z, xt, xs0, xs1, ts0, ts1)*tpss_aux(z, xt, ts0, ts1)^3):
