(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

alpha := (4/(9*Pi))^(1/3):

# Eq. (46) of Phys. Rev. A 73, 032506 (2006)
# I am not sure if this is the on-top pair correlation functional
# intended by the authors
g0_C    :=  0.08193:
g0_D    := -0.01277:
g0_E    :=  0.001859:
g0_dd   :=  0.7524:
g0_aHD  := -alpha*(Pi^2 + 6*log(2) - 3)/(5 * Pi):
g0_B    := -2*g0_aHD - g0_dd:
pmgb_g0 := rs ->
        (1 - g0_B*rs + g0_C*rs^2 + g0_D*rs^3 + g0_E*rs^4)*exp(-g0_dd*rs)/2:

pmgb_gc0 := rs -> pmgb_g0(rs) - 1/2:

# Eq. (14)
pmgb_phi := (n, z) ->
         (opz_pow_n(z,n/3) + opz_pow_n(-z,n/3))/2:

# Eq. (22)
pmgb_Q_a := 5.84605:
pmgb_Q_c := 3.91744:
pmgb_Q_d := 3.44851:
pmgb_Q_b := pmgb_Q_d - 3*Pi*alpha/(4*log(2) - 4):

pmgb_Q := x ->
       ((2*log(2) - 2)/Pi^2) * log(
        (1 + pmgb_Q_a*x + pmgb_Q_b*x^2 + pmgb_Q_c*x^3)/
        (1 + pmgb_Q_a*x + pmgb_Q_d*x^2)
       ):

# Eq. (31)
pmgb_b0 := (rs) -> 0.784949*rs:

# Eq. (32)
pmgb_gpp1 := (rs) -> (2^(5/3)/(5*alpha^2*rs^2)) *
          (1 - 0.02267*rs)/(1 + 0.4319*rs + 0.04*rs^2):

# Eq. (33)
pmgb_D2 := (rs) -> (-0.388*rs + 0.676*rs^2)*exp(-0.547*rs)/rs^2:

# Eq. (34)
pmgb_D3 := (rs) -> (-4.95*rs + rs^2)*exp(-0.31*rs)/rs^3:

# Eq. (28)
pmgb_cc4 := (rs, z) ->
       + (opz_pow_n( z,2)/4) * pmgb_gpp1(rs*(2/(1 + z))^(1/3))
       + (opz_pow_n(-z,2)/4) * pmgb_gpp1(rs*(2/(1 - z))^(1/3))
       + (1 - z^2) * pmgb_D2(rs)
       - pmgb_phi(8, z)/(5*alpha^2*rs^2):

# Eq. (29)
pmgb_cc5 := (rs, z) ->
       + (opz_pow_n( z,2)/4) * pmgb_gpp1(rs*(2/(1 + z))^(1/3))
       + (opz_pow_n(-z,2)/4) * pmgb_gpp1(rs*(2/(1 - z))^(1/3))
       + (1 - z^2) * pmgb_D3(rs):

# Eq. (30)
pmgb_C2 := (rs, z) ->
        -3*(1 - z^2)*pmgb_gc0(rs)/(8*rs^3):

pmgb_C3 := (rs, z) ->
        -(1 - z^2)*pmgb_g0(rs)/(sqrt(2*Pi)*rs^3):

pmgb_C4 := (rs, z) ->
        -9*pmgb_cc4(rs, z)/(64*rs^3):

pmgb_C5 := (rs, z) ->
        -9*pmgb_cc5(rs, z)/(40*sqrt(2*Pi)*rs^3):

pmgb_a1 := (rs, z) ->
        4*pmgb_b0(rs)^6*pmgb_C3(rs, z) + pmgb_b0(rs)^8*pmgb_C5(rs, z):

pmgb_a2 := (rs, z) ->
        4*pmgb_b0(rs)^6*pmgb_C2(rs, z) + pmgb_b0(rs)^8*pmgb_C4(rs, z)
        + 6*pmgb_b0(rs)^4*f_pw(rs, z):

pmgb_a3 := (rs, z) ->
        pmgb_b0(rs)^8*pmgb_C3(rs, z):

pmgb_a4 := (rs, z) ->
        pmgb_b0(rs)^8*pmgb_C2(rs, z) + 4*pmgb_b0(rs)^6*f_pw(rs, z):

pmgb_a5 := (rs, z) ->
        pmgb_b0(rs)^8*f_pw(rs, z):

# Eq. (26)
pmgb_ec_LR := (rs, z) -> (
           + pmgb_phi(2, z)^3*pmgb_Q(p_a_hyb_omega_0_*sqrt(rs)/pmgb_phi(2, z))
           + pmgb_a1(rs, z)*p_a_hyb_omega_0_^3
           + pmgb_a2(rs, z)*p_a_hyb_omega_0_^4
           + pmgb_a3(rs, z)*p_a_hyb_omega_0_^5
           + pmgb_a4(rs, z)*p_a_hyb_omega_0_^6
           + pmgb_a5(rs, z)*p_a_hyb_omega_0_^8
           )/(1 + pmgb_b0(rs)^2*p_a_hyb_omega_0_^2)^4:

f := (rs, z) -> pmgb_ec_LR(rs, z):
