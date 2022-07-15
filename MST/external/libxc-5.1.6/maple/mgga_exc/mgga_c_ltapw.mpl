(* type: mgga_exc *)
(* prefix:
  mgga_c_ltapw_params *params;

  assert(p->params != NULL);
  params = (mgga_c_ltapw_params * )(p->params);
*)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

(* kinetic energy density to electron density *)
nt_tau := t -> (t/K_FACTOR_C)^(3*params_a_ltafrac/5):

(* effective density *)
n_eff_s := (rs, z, t) -> n_spin(rs, z) * nt_tau(t):
n_eff   := (rs, z, ts0, ts1) -> n_eff_s(rs, z, ts0) + n_eff_s(rs, -z, ts1):

(* recompute rs and zeta *)
eff_rs := (rs, z, ts0, ts1) -> r_ws(n_eff(rs, z, ts0, ts1)):
eff_z  := (rs, z, ts0, ts1) ->
  (n_eff_s(rs, z, ts0) - n_eff_s(rs, -z, ts1))/n_eff(rs, z, ts0, ts1):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  f_pw(eff_rs(rs, z, ts0, ts1), eff_z(rs, z, ts0, ts1)):

