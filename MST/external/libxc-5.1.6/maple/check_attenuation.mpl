(* 2020 Susi Lehtola and Miguel A. L. Marques

   This script is used to check the convergence of the series
   expansions of the various attenuation functions, and to
   establish cutoff values for their series expansions.
*)

$include "attenuation.mpl"

(* For exact results we use 1000 digit *)
exact_digits := 1000:
(* Double precision has 15 digits *)
double_digits := 15:
(* An acceptable relative error at double precision is 1e-13 *)
error_thr := 1e-13:

check_asymptotics := proc(f, a);
    (* Grid spacing *)
    da := 0.01:
    (* Value for a to start with *)
    amax := 5:

    (* First, we compare we find a point where to make the cutoff *)
    for acut from amax by -da to da do
        (* Exact result *)
        Digits := exact_digits:
        exact := evalf(f(acut)):
        (* Double precision *)
        Digits := double_digits:
        local doubleprec := evalf(f(acut)):
        (* Error *)
        local err := doubleprec/exact-1:
        #printf("Cutoff %e exact % e double % e error % e\n", acut, exact, doubleprec, err):
        if (abs(err) < error_thr) then
            break:
        end if:
    end do:

    (* Now we find the series expansion that has the same level of agreement *)
    for expord from 4 to 1000 by 2 do
        f_series := a -> eval(throw_out_large_n(convert(series(f(b), b=infinity, expord+padding_order), polynom), expord), b=a):
        local ser := evalf(f_series(acut)):
        local err := ser/exact-1:
        #printf("Expansion order %3d original % e asymptotic % e error % "
        #       "e\n", expord, exact, ser, err);
        if (abs(err) < error_thr) then
            (* Check if the expansion is accurate everywhere *)
            accurate := true:
            for aval from acut by da to amax do
                (* Exact result *)
                Digits := exact_digits:
                lexact := evalf(f(aval)):
                (* Double precision *)
                Digits := double_digits:
                lser := evalf(f_series(aval)):
                (* Error *)
                local lerr := lser/lexact-1:
                #printf("a= %e exact % e double % e error % e\n", aval, lexact, lser, lerr):
                if (abs(lerr) > error_thr) then
                    accurate := false:
                    break:
                end if:
            end do:
            if accurate then
                break:
            end if:
        end if:
    end do:

    printf("Cutoff should be at a= %.2f and use a series expansion of "
           "order %d\n", acut, expord):
end proc:

printf("attenuation_erf\n"):
check_asymptotics(attenuation_erf0, a):

printf("\nattenuation_erf_f2\n"):
check_asymptotics(attenuation_erf_f20, a):

printf("\nattenuation_erf_f3\n"):
check_asymptotics(attenuation_erf_f30, a):

printf("\nattenuation_gau\n"):
check_asymptotics(attenuation_gau0, a):

printf("\nattenuation_yukawa\n"):
check_asymptotics(attenuation_yukawa0, a):

