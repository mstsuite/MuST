(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

(* This was checked against code kindly provided by Emil Proynov *)

(* Equations (24) *)
a_24_i := [-113.693369789727190, 24.00502151278711440, 49.34131295839670750,
            -23.8242372168379302, 0.944080741695104794, 0.000293039144178338]:
b_24_i := [-109.74263493216910, 16.2663129444242415, 54.4034331373908366,
            -25.154009904187990, 1.0]:

f_r := rs -> add(a_24_i[i]*rs^(i-1), i=1..6)/add(b_24_i[i]*rs^(i-1), i=1..5):

(* Equations (25) *)
c_25_i := [-0.32481568604919886, 1.180131465463191050, -1.42693041498421640,
            0.580344063812247980, -0.01099122367291440]:
d_25_i := [-0.57786103193239430, 2.09708505883490736, -2.52188183586948180, 1.0]:

f_s := z -> add(c_25_i[i]*z^(i-1), i=1..5)/add(d_25_i[i]*z^(i-1), i=1..4):

(* Equation (23) *)
(* The factor 1.28 is absent from the paper, but it is in the original code. See erratum *)
ss := (rs, z) -> f_r(rs)*f_s(z)*1.28:

(* Equation (22) *)
alpha_z := (rs, z) -> 2/(opz_pow_n(z,ss(rs, z)) + opz_pow_n(-z,ss(rs, z))):

(* Equation (21) *)
eta6  := 0.41081146652128:
eta7  := 0.599343256903515:
eta8  := 1.70939476802168:
eta9  := 0.077123208419481:
eta10 := 0.46958449007619:

alpha_n := rs ->
  + eta6
  + eta7*exp( -eta8*rs^(1/3))*rs^(2/3)
  + eta9*exp(-eta10*rs^(1/3))*rs^(1/3):

(* Equation (20) *)
alpha_eff := (rs, z) -> alpha_n(rs)*alpha_z(rs, z):

(* Equation (19) *)
eta1 :=  0.538074483500437:
eta2 := -2.226094990985190:
eta3 :=  0.837303782322808:
eta4 :=  2.619709858963178:
eta5 :=  1.036657594643520:

beta_eff := rs ->
  + eta1
  + eta2*exp(-eta3*rs^(1/3))*rs^(1/4)
  + eta4*exp(-eta5*rs^(1/3))*rs^(1/3):

(* Equation (15), see erratum *)
ax   := (3*Pi^2)^(1/3):
k_fs := (rs, z) -> ax*RS_FACTOR/rs * opz_pow_n(z,1/3):

(* Equation (17) *)
k_uu := (rs, z) -> alpha_eff(rs,  z)*k_fs(rs,  z):
k_dd := (rs, z) -> alpha_eff(rs, -z)*k_fs(rs, -z):

(* Equation (18) *)
k_ud := (rs, z) -> beta_eff(rs)
  * 2*k_fs(rs, z)*k_fs(rs, -z)/(k_fs(rs, z) + k_fs(rs, -z)):

(* Table III *)
a1  := 0.1846304394851914:
a2  := 5.93965654951900799:
a3  := 2.36958012866641818:
a4  := .51188865525958770e-1:
a5  := .9576892532004281e-1:
a6  := .283592616144882565e-1:
a7  := .226274169979695208e-1:
a8  := .531736155271654809e-2:
a9  := .1915378506400854:
a10 := .1473137771194929:
a11 := .1528250938350897:
a12 := 1.01508307543839117:
a13 := .7641254691754473e-1:
a14 := .898537460263473410:
a15 := .1795667349750801e-1:
a16 := .3461820740347690e-1:
a17 := .3591334699501599e-1:
a18 := .222017353476155799:

c1  := 132.479090287794355:
c2  := 32.4014708516771368:
c3  := 22.5664453162503806:
c4  := 11.2832226581251903:
c5  := .401060523940960082:
c6  := evalf(0.32):
c7  := .751988482389300153e-1:
c8  := 116.935042647480910:
c9  := 29.6240023046901289:
c10 := .482257181994472723:
c11 := .246903981179097557:
c12 := evalf(1/2):
c13 := .410709696778185459:
c14 := .105323524476768857:
c15 := 14.5650971711659670:
c16 := .781250000000000000:
c17 := .623347313127238558:
c18 := .146484374999999999:
c19 := 111.811548105797788:
c20 := .160041105570901272:
c21 := evalf(.78125):
c22 := .32086695060795739:
c23 := 13.2844495072998436:
c24 := .268418671319107341:
c25 := .471060597934991862:
c26 := evalf(1/4):
c27 := .252882919616989509:
c28 := .720485831127149779e-1:
c29 := 42.6490544891031073:

(* Definitions in the beginning of the appendix *)
D_1 := k -> a6*k^2 + a7*k + a8:
D_2 := k -> a1*k^2 + a10*k + a16:
D_3 := k -> a5*k^2 + a13*k + a15:
D_4 := k -> a9*k^2 + a11*k + a17:
D_5 := k -> c5*k^2 + c6*k + c7:
D_6 := k -> c12*k^2 + c13*k + c14:
D_7 := k -> c16*k^2 + c17*k + c18:
D_8 := k -> sqrt(c26*k^2 + c27*k + c28):

(* Equation (10) *)
Q_1ud := k ->  1/D_1(k) * (
  - arctan(a2*k + a3)*D_2(k)/k - log(D_1(k))*D_3(k)/k
  + log(k)*D_4(k)/k - a4*k + a12 + a14/k + a18/k^2
):

(* Equation (11) *)
Q_2ud := k ->
  - c1/k - c2/k^2 - c3*log(k)/k + c4*log(D_5(k))/k
  + c8*arctan(a2*k + a3)/k + c9*log(k + c10)/k - c11/k*log(D_6(k)):

(* Equation (12) *)
Q_3ud := k ->
  + c19*arctan(c20/(c21*k + c22))/k - c23*arctanh((c24 + c25*k)/D_8(k))/k
  - c15*log(D_7(k))/k - c29*D_8(k)/k^2:

(* Equation (9) *)
ec_opp := (rs, z) ->
  (1 - z^2)/4*(Q_1ud(k_ud(rs, z)) + Q_2ud(k_ud(rs, z)) + Q_3ud(k_ud(rs, z))):

(* Equation (13) *)
ec_par := (rs, z) ->
  + opz_pow_n( z,2)/8*(Q_1ud(k_uu(rs, z)) + Q_2ud(k_uu(rs, z)) + Q_3ud(k_uu(rs, z)))
  + opz_pow_n(-z,2)/8*(Q_1ud(k_dd(rs, z)) + Q_2ud(k_dd(rs, z)) + Q_3ud(k_dd(rs, z))):

f := (rs, z) -> n_total(rs)*(ec_opp(rs, z) + ec_par(rs, z)):
