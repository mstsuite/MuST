(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

# Minimum and maximum functions that I can differentiate

m_min := (x1, x2) -> my_piecewise3(x1 > x2, x2, x1):
m_max := (x1, x2) -> my_piecewise3(x1 > x2, x1, x2):
m_abs := (x)      -> my_piecewise3(x > 0, x, -x):

# Teach maple to differentiate piecewise functions

`diff/my_piecewise3` :=
    proc(c, x1, x2, x) my_piecewise3(c, diff(x1, x), diff(x2, x)) end proc:

`diff/my_piecewise5` :=
    proc(c1, x1, c2, x2, x3, x) my_piecewise5(c1, diff(x1, x), c2, diff(x2, x), diff(x3,x)) end proc:

# This is the derivative of xc_E1_scaled = -exp(x)*Ei(-x) = exp(x)E1(x)

`diff/xc_E1_scaled` :=
  proc(y, x) (xc_E1_scaled(y) - 1/y) * diff(y, x) end proc:

# The derivative of xc_erfcx = exp(x^2)*erfc(x) = 2*x*xc_erfcx(x) - 2/sqrt(Pi)
`diff/xc_erfcx` :=
  proc(y, x) (2*y*xc_erfcx(y) - 2/sqrt(Pi)) * diff(y, x) end proc:


# a series of useful definitions

M_C         := 137.0359996287515: (* speed of light *)

X2S         := 1/(2*(6*Pi^2)^(1/3)):
X2S_2D      := 1/(2*(4*Pi)^(1/2)):

X_FACTOR_C    := 3/8*(3/Pi)^(1/3)*4^(2/3):
X_FACTOR_2D_C := 8/(3*sqrt(Pi)):
K_FACTOR_C    := 3/10*(6*Pi^2)^(2/3):

MU_GE       := 10/81:
MU_PBE      := 0.06672455060314922*(Pi^2)/3:
KAPPA_PBE   := 0.8040:

# generic conversion functions

$ifdef xc_dimensions_1d
DIMENSIONS   := 1:
RS_FACTOR    := 1/2:
$elif xc_dimensions_2d
DIMENSIONS   := 2:
RS_FACTOR    := 1/sqrt(Pi):
LDA_X_FACTOR := -X_FACTOR_2D_C:
$else
DIMENSIONS   := 3:
RS_FACTOR    := (3/(4*Pi))^(1/3):
LDA_X_FACTOR := -X_FACTOR_C:
$endif

r_ws       := n  -> RS_FACTOR/n^(1/DIMENSIONS):
n_total    := rs -> (RS_FACTOR/rs)^DIMENSIONS:

n_spin     := (rs, z) -> simplify((1 + z)*n_total(rs)/2):
sigma_spin := (rs, z, xs) -> simplify(xs^2*n_spin(rs, z)^(8/3)):
t_total    := (z, ts0, ts1) ->
  (ts0*((1 + z)/2)^(5/3) + ts1*((1 - z)/2)^(5/3)):
u_total    := (z, us0, us1) -> t_total(z, us0, us1):

# useful formulas that enter several functionals follow

# von Weizsäcker term
t_vw := (z, xt, us0, us1) -> (xt^2 - u_total(z, us0, us1))/8:

# Screening for extreme values of zeta
z_thr := z -> my_piecewise5(
  1 + z <= p_a_zeta_threshold, p_a_zeta_threshold - 1,
  1 - z <= p_a_zeta_threshold, 1 - p_a_zeta_threshold,
  z):

opz_pow_n := (z, n) -> my_piecewise3(1 + z <= p_a_zeta_threshold, (p_a_zeta_threshold)^n, (1 + z)^n):

# See Eq. (9) of Perdew1992_13244
f_zeta    := z -> (opz_pow_n(z,4/3) + opz_pow_n(-z,4/3) - 2)/(2^(4/3) - 2):
f_zeta_2d := z -> 1/2*(opz_pow_n(z,3/2) + opz_pow_n(-z,3/2)):

# used in several correlation functionals
mphi := z -> (opz_pow_n(z,2/3) + opz_pow_n(-z,2/3))/2:
tt   := (rs, z, xt) -> xt/(4*2^(1/3)*mphi(z)*sqrt(rs)):

# in the paper it is beta_a = 0.066725
beta_Hu_Langreth := rs -> 0.066724550603149220*(1 + 0.1*rs)/(1 + 0.1778*rs):

# Generate exchange and kinetic functionals from the expression for the
# enhancement factor
lda_x_spin := (rs, z) -> simplify(LDA_X_FACTOR*opz_pow_n(z,1 + 1/DIMENSIONS)*2^(-1-1/DIMENSIONS)*(RS_FACTOR/rs)):
lda_k_spin := (rs, z) -> simplify(K_FACTOR_C*opz_pow_n(z,5/3)*2^(-5/3)*(RS_FACTOR/rs)^2):

# Screen out small densities as well as the zero component from fully spin polarized densities due to terms of the form (1+z)^{-n}
screen_dens := (rs, z) -> (n_spin(rs,  z) <= p_a_dens_threshold):

# For most functionals zeta screening occurs inside the funcitonal,
#  but for B88 and B94 correlation need to screen out also outside
screen_dens_zeta := (rs, z) -> screen_dens(rs, z) or (1 + z <= p_a_zeta_threshold):

# non-separable GGA exchange
gga_exchange_nsp := (func, rs, z, xs0, xs1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, lda_x_spin(rs, z_thr( z))*func(rs, z_thr( z), xs0))
  + my_piecewise3(screen_dens(rs, -z), 0, lda_x_spin(rs, z_thr(-z))*func(rs, z_thr(-z), xs1)):
# GGA exchange
gga_exchange := (func, rs, z, xs0, xs1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, lda_x_spin(rs, z_thr( z))*func(xs0))
  + my_piecewise3(screen_dens(rs, -z), 0, lda_x_spin(rs, z_thr(-z))*func(xs1)):
# GGA kinetic energy
gga_kinetic := (func, rs, z, xs0, xs1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, lda_k_spin(rs, z_thr( z))*func(xs0))
  + my_piecewise3(screen_dens(rs, -z), 0, lda_k_spin(rs, z_thr(-z))*func(xs1)):

# non-separable meta-GGA exchange
mgga_exchange_nsp := (func, rs, z, xs0, xs1, u0, u1, t0, t1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, lda_x_spin(rs, z_thr( z))*func(rs, z_thr( z), xs0, u0, t0))
  + my_piecewise3(screen_dens(rs, -z), 0, lda_x_spin(rs, z_thr(-z))*func(rs, z_thr(-z), xs1, u1, t1)):
# meta-GGA exchange
mgga_exchange := (func, rs, z, xs0, xs1, u0, u1, t0, t1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, lda_x_spin(rs, z_thr( z))*func(xs0, u0, t0))
  + my_piecewise3(screen_dens(rs, -z), 0, lda_x_spin(rs, z_thr(-z))*func(xs1, u1, t1)):
# meta-GGA kinetic energy
mgga_kinetic := (func, rs, z, xs0, xs1, u0, u1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, lda_k_spin(rs, z_thr( z))*func(xs0, u0))
  + my_piecewise3(screen_dens(rs, -z), 0, lda_k_spin(rs, z_thr(-z))*func(xs1, u1)):

# This is the Stoll decomposition in our language
lda_stoll_par  := (lda_func, rs, z) ->
  my_piecewise3(screen_dens_zeta(rs,  z), 0, opz_pow_n(z, 1)/2 * lda_func(rs*2^(1/3)*opz_pow_n(z, -1/3), 1)):

lda_stoll_perp := (lda_func, rs, z) ->
  + lda_func(rs, z)
  - lda_stoll_par(lda_func, rs,  z)
  - lda_stoll_par(lda_func, rs, -z):

gga_stoll_par  := (gga_func, rs, z, xs, spin) ->
  my_piecewise3(screen_dens_zeta(rs, z), 0,
    gga_func(rs*2^(1/3)*opz_pow_n(z,-1/3), spin, xs, xs*(1 + spin)/2, xs*(1 - spin)/2) * opz_pow_n(z,1)/2):

# Curvature of the Fermi hole (without the current term)
Fermi_D := (xs, ts) -> 1 - xs^2/(8*ts):

# correction to Fermi_D similar to the one found in
#  JCP 127, 214103 (2007); doi: http://dx.doi.org/10.1063/1.2800011
Fermi_D_corrected := (xs, ts) -> (1 - xs^2/(8*ts)) * (1 - exp(-4*ts^2/params_a_Fermi_D_cnst^2)):

# Becke function used in several correlation functionals
b88_R_F := (f_x, rs, z, xs) -> 1/(2*X_FACTOR_C*n_spin(rs, z)^(1/3)*f_x(xs)):
b88_zss := (css, f_x, rs, z, xs) -> 2*css*b88_R_F(f_x, rs, z, xs):
b88_zab := (cab, f_x, rs, z, xs0, xs1) -> cab*(
  + my_piecewise3(screen_dens(rs,  z), 0, b88_R_F(f_x, rs, z_thr( z), xs0))
  + my_piecewise3(screen_dens(rs, -z), 0, b88_R_F(f_x, rs, z_thr(-z), xs1))
  ):

# The meta-GGA version
b94_R_F := (f_x, rs, z, xs, us, ts) -> 1/(2*X_FACTOR_C*n_spin(rs, z)^(1/3)*f_x(xs,us,ts)):
b94_zss := (css, f_x, rs, z, xs, us, ts) -> 2*css*b94_R_F(f_x, rs, z, xs, us, ts):
b94_zab := (cab, f_x, rs, z, xs0, xs1, us0, us1, ts0, ts1) -> cab*(
  + my_piecewise3(screen_dens(rs,  z), 0, b94_R_F(f_x, rs, z_thr( z), xs0, us0, ts0))
  + my_piecewise3(screen_dens(rs, -z), 0, b94_R_F(f_x, rs, z_thr(-z), xs1, us1, ts1))
  ):

# Power series often used in mggas
mgga_w := t -> (K_FACTOR_C - t)/(K_FACTOR_C + t):
mgga_series_w := (a, n, t) -> add(a[i]*mgga_w(t)^(i-1), i=1..n):

# Used in screened functionals
kF := (rs, z) -> (3*Pi^2)^(1/3) * opz_pow_n(z,1/3) * RS_FACTOR/rs:
nu := (rs, z) -> p_a_hyb_omega_0_/kF(rs, z):


(* Maple 2020 has a bug where series aren't computed to the order requested; this circumvents that *)
padding_order := 30:

(* Cap polynomial expansion to given order (throws out the padded terms, if any) *)
throw_out_large_n := proc(X,n) select(t -> abs(degree(t,{b}))<=n, X); end proc:

(* Function that makes f(a) smooth for a->infty *)
enforce_smooth_lr := proc(f, a, a_cutoff, expansion_order);
  (* Calculate large-a expansion *)
  f_large := a -> eval(throw_out_large_n(convert(series(f(b), b=infinity, expansion_order+padding_order), polynom), expansion_order), b=a):
  (* Return the series expansion for large a; also remove any numerical overflows from the original branch  *)
  my_piecewise3(a >= a_cutoff, f_large(m_max(a, a_cutoff)), f(m_min(a, a_cutoff))):
end proc:
