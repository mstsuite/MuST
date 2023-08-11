(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$include "scaling.mpl"

(* parameters *)
wpbeh_A :=  1.0161144:
wpbeh_B := -0.37170836:
wpbeh_C := -0.077215461:
wpbeh_D :=  0.57786348:
wpbeh_E := -0.051955731:

(*
 Note that kF has a 6 and not a 3 as it should in principle
 be. This is because the HSE formula, if one would take the papers
 seriously, does not fulfill the spin sum-rule. This is probably
 an oversight from them. So, we have to choose, either a 6 or a 3.

 Nwchem seems to have the factor of 6, but VASP and espresso have
 a 3. This would amount to rescaling omega by a factor of
 cbrt(2). We follow the quantum chemistry community and put the 6.
*)

(* Cutoff criterion below which to use polynomial expansion *)
EGscut     := 0.08:
wcutoff    := 14:
expfcutoff := 700:

(* first let us calculate H(s) *)
wpbeh_Ha1 := 0.00979681:
wpbeh_Ha2 := 0.0410834:
wpbeh_Ha3 := 0.187440:
wpbeh_Ha4 := 0.00120824:
wpbeh_Ha5 := 0.0347188:
wpbeh_H := s ->
  + (wpbeh_Ha1*s^2 + wpbeh_Ha2*s^4)
  / (1 + wpbeh_Ha3*s^4 + wpbeh_Ha4*s^5 + wpbeh_Ha5*s^6):

(*
  Now we calculate F(s). We use the parameters that were in the original code,
  but these constants are:

  wpbeh_Fc1 := 4*wpbeh_A^2/(9*wpbeh_C) + (wpbeh_B - wpbeh_A*wpbeh_D)/wpbeh_C:
  wpbeh_Fc2 := -4/(3*36*wpbeh_C):
*)
wpbeh_Fc1 := 6.4753871:
wpbeh_Fc2 := 0.47965830:
wpbeh_F := s ->
  wpbeh_Fc1*wpbeh_H(s) + wpbeh_Fc2:

(* several auxiliary variables *)
eb1  := w -> my_piecewise3(w < wcutoff, 1.455915450052607, 2):
aux1 := s -> wpbeh_D + s^2*wpbeh_H(s):
aux2 := s -> 9*wpbeh_H(s)*s^2/(4*wpbeh_A):
aux3 := (w, s) -> aux1(s) + w^2:
aux4 := (w, s) -> s^2*wpbeh_H(s) + eb1(w)*w^2:
aux5 := (w, s) -> 9*aux4(w, s)/(4*wpbeh_A):
aux6 := (w, s) -> wpbeh_D + aux4(w, s):

(* and now G(s) *)
Ga := s ->
  + sqrt(Pi)*(
    + 15*wpbeh_E
    + 6*wpbeh_C*(1 + wpbeh_F(s)*s^2)*aux1(s)
    + 4*wpbeh_B*aux1(s)^2
    + 8*wpbeh_A*aux1(s)^3)
  / (16*aux1(s)^(7/2))
  - (3*Pi/4)*sqrt(wpbeh_A)*exp(aux2(s))*(1 - erf(sqrt(aux2(s)))):
Gb := s ->
  15*sqrt(Pi)*s^2/(16*aux1(s)^(7/2)):

wpbeh_EGa1 := -0.02628417880:
wpbeh_EGa2 := -0.07117647788:
wpbeh_EGa3 :=  0.08534541323:
wpbeh_EG := s -> my_piecewise3(s > EGscut,
  -(3*Pi/4 + Ga(s))/Gb(s),
  wpbeh_EGa1 + wpbeh_EGa2*s^2 + wpbeh_EGa3*s^4
):

term2 := s-> (
  + aux1(s)^2*wpbeh_B
  + aux1(s)*wpbeh_C
  + 2*wpbeh_E
  + aux1(s)*s^2*wpbeh_C*wpbeh_F(s)
  + 2*s^2*wpbeh_EG(s)
)/(2*aux1(s)^3):

term3 := (w, s) -> -w*(
  + 4*aux3(w, s)^2*wpbeh_B
  + 6*aux3(w, s)*wpbeh_C
  + 15*wpbeh_E
  + 6*aux3(w, s)*s^2*wpbeh_C*wpbeh_F(s)
  + 15*s^2*wpbeh_EG(s)
)/(8*aux1(s)*aux3(w, s)^(5/2)):

term4 := (w, s) -> -w^3*(
  + aux3(w, s)*wpbeh_C
  + 5*wpbeh_E
  + aux3(w, s)*s^2*wpbeh_C*wpbeh_F(s)
  + 5*s^2*wpbeh_EG(s)
)/(2*aux1(s)^2*aux3(w, s)^(5/2)):

term5 := (w, s) -> -w^5*(
  + wpbeh_E
  + s^2*wpbeh_EG(s)
)/(aux1(s)^3*aux3(w, s)^(5/2)):

t10 := (w, s) ->
  1/2*wpbeh_A*log(aux4(w, s)/aux6(w, s)):

(* Use simple gaussian approximation for large w *)
term1_largew := (w, s) ->
  -1/2*wpbeh_A*(-xc_E1_scaled(aux5(w, s)) + log(aux6(w, s)) - log(aux4(w, s))):

(* For everything else use the full blown expression *)
ea1 := -1.128223946706117:
ea2 :=  1.452736265762971:
ea3 := -1.243162299390327:
ea4 :=  0.971824836115601:
ea5 := -0.568861079687373:
ea6 :=  0.246880514820192:
ea7 := -0.065032363850763:
ea8 :=  0.008401793031216:

np1 := w ->
  - 1.5*ea1*sqrt(wpbeh_A)*w
  + 27*ea3*w^3/(8*sqrt(wpbeh_A))
  - 243*ea5*w^5/(32*(wpbeh_A)^(3/2))
  + 2187*ea7*w^7/(128*(wpbeh_A)^(5/2)):

np2 := w ->
  - wpbeh_A
  + 9*ea2*w^2/4.0
  - 81*ea4*w^4/(16*wpbeh_A)
  + 729*ea6*w^6/(64*wpbeh_A^2)
  - 6561*ea8*w^8/(256*wpbeh_A^3):

t1 := (w, s) ->
  1/2*(np1(w)*Pi*xc_erfcx(sqrt(aux5(w, s))) - np2(w)*xc_E1_scaled(aux5(w, s))):

f2 := (w, s) ->
  1/2*ea1*sqrt(Pi)*wpbeh_A/sqrt(aux6(w, s)):
f3 := (w, s) ->
  1/2*ea2*wpbeh_A/aux6(w, s):
f4 := (w, s) ->
  ea3*sqrt(Pi)*(-9/(8*sqrt(aux4(w, s))) + 0.25*wpbeh_A/aux6(w, s)^(3/2)):
f5 := (w, s) ->
  (ea4/128)*(-144/aux4(w, s) + 64*wpbeh_A/aux6(w, s)^2):
f6 := (w, s) ->
  ea5*(3*sqrt(Pi)*(3*aux6(w, s)^(5/2)*(9*aux4(w, s) - 2*wpbeh_A)
    + 4*aux4(w, s)^(3/2)*wpbeh_A^2))/(32*aux6(w, s)^(5/2)*aux4(w, s)^(3/2)*wpbeh_A):
f7 := (w, s) ->
  ea6*((32*wpbeh_A/aux6(w, s)^3 + (-36 + 81*s^2*wpbeh_H(s)/wpbeh_A)/aux4(w, s)^2))/32:
f8 := (w, s) ->
  ea7*(-3*sqrt(Pi)*(-40*aux4(w, s)^(5/2)*wpbeh_A^3
    + 9*aux6(w, s)^(7/2)*(27*aux4(w, s)^2 - 6*aux4(w, s)*wpbeh_A + 4*wpbeh_A^2)))
  /(128*aux6(w, s)^(7/2)*aux4(w, s)^(5/2)*wpbeh_A^2):
f9 := (w, s) -> (
  + 324*ea6*eb1(w)*aux6(w, s)^4*aux4(w, s)*wpbeh_A
  + ea8*(384*aux4(w, s)^3*wpbeh_A^3
    + aux6(w, s)^4*(-729*aux4(w, s)^2 + 324*aux4(w, s)*wpbeh_A - 288*wpbeh_A^2))
  )/(128*aux6(w, s)^4*aux4(w, s)^3*wpbeh_A^2):

t2t9 := (w, s) ->
  + f2(w, s)*w + f3(w, s)*w^2 + f4(w, s)*w^3 + f5(w, s)*w^4
  + f6(w, s)*w^5 + f7(w, s)*w^6 + f8(w, s)*w^7 + f9(w, s)*w^8:

term1 := (w, s) -> my_piecewise3(
  w > wcutoff, term1_largew(w, s),
  t1(m_min(w, wcutoff), s) + t2t9(m_min(w, wcutoff), s) + t10(m_min(w, wcutoff), s)
):

f_wpbeh0 := (w, s) -> - 8/9 *(
  term1(w, s) + term2(s) + term3(w, s) + term4(w, s) + term5(w, s)
):

f_wpbeh := (rs, z, x) ->
  f_wpbeh0(nu(rs, z), m_max(1e-15, s_scaling_2(X2S*x))):

f  := (rs, z, xt, xs0, xs1) ->
  gga_exchange_nsp(f_wpbeh, rs, z, xs0, xs1):
