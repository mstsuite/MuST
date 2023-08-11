(* type: lda_exc *)

a :=  0.0311:
b := -0.048:
c :=  0.009:
d := -0.017:

f_ls := (rs, zeta) -> a*log(rs) + b + c*rs*log(rs) + d*rs:
f_ms := (rs, zeta) -> -0.06156 + 0.01898*log(rs):
f_hs := (rs, zeta) -> 0.438/rs + 1.325/rs^(3/2) - 1.47/rs^2 - 0.4/rs^(5/2):

f := (rs, zeta) -> my_piecewise3(
  rs < 0.7, f_ls(rs, zeta),  my_piecewise3(rs < 10,  f_ms(rs, zeta), f_hs(rs, zeta))):

