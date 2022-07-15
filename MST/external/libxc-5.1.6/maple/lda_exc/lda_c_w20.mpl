(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

(* Data in table 1 *)
w20_a0 := proc(z)
  if (z=0) then
    (1-log(2))/Pi^2
  elif (z=1) then
    (1-log(2))/(2*Pi^2)
  end if
end proc:

w20_a1 := proc(z)
  if (z=0) then
    (9*Pi/4)^(-1/3) * 1/(4*Pi^3) * ( 7*Pi^2/6 - 12*log(2) - 1 )
  elif (z=1) then
    2^(-4/3) * (9*Pi/4)^(-1/3) * 1/(4*Pi^3) * ( 13*Pi^2/12 - 12*log(2) + 1/2 )
  end if
end proc:

w20_b0 := proc(z)
  if (z=0) then
    -0.071100 + log(2)/6 - 3/(4*Pi^2)*evalf(Zeta(3))
  elif (z=1) then
    -0.049917 + log(2)/6 - 3/(4*Pi^2)*evalf(Zeta(3))
  end if
end proc:

w20_b1 := proc(z)
  if (z=0) then
    -0.01
  elif (z=1) then
    (* Compare eq 12 to eqs 15 and 16 *)
    0
  end if
end proc:

w20_f0 := -0.9:
w20_f1 := 1.5:
w20_f2 := 0:

(* eq 6 *)
w20_cs := z -> 3/10*(9*Pi/4)^(2/3) * 1/2*((1+z)^(5/3) + (1-z)^(5/3)):
(* eq 7 *)
w20_cx := z -> -3/(4*Pi)*(9*Pi/4)^(1/3) * 1/2*((1+z)^(4/3) + (1-z)^(4/3)):

(* eq 8 *)
w20_ec := (rs, z) -> -w20_a0(z)/2 * log(1 + w20_DF(rs,z,w20_f0-w20_cx(z))/rs + w20_E(rs,z)/rs^(3/2) + w20_DF(rs,z,w20_f2-w20_cs(z))/rs^2) + w20_G(rs,z):

(* eqs 9 and 11 only differ by the f_i - c_j(z) term *)
w20_DF := (rs, z, cfterm) -> exp(-2*w20_b0(z)/w20_a0(z)) - 2 * (1 - exp(-(rs/100)^2))*( cfterm/w20_a0(z) + 1/2 * exp(-2*w20_b0(z)/w20_a0(z)) ):

(* eq 10 *)
w20_E := (rs, z) -> - 2*(1 - exp(-(rs/100)^2))*w20_f1 / w20_a0(z):

(* eq 12, rewritten in terms of a decaying exponential to avoid overflow *)
w20_G := (rs, z) -> rs*exp(-(rs/100)^2) / (exp(-(rs/100)^2) + 10*rs^(5/4)) * ( -w20_a1(z)*log(1 + 1/rs) + w20_b1(z) ):

(* eq 17 *)
f_w20 := (rs, zeta) -> w20_ec(rs,0) + (w20_ec(rs,1) - w20_ec(rs,0))*f_zeta(zeta):

f := (rs, zeta) -> f_w20(rs, zeta):
