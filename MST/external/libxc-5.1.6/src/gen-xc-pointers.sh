#!/bin/bash
#
# Copyright (C) 2020 Susi Lehtola
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

ingr=(rho sigma lapl tau)
der=("VXC" "FXC" "KXC" "LXC")

max=4

for((mode=0;mode<2;mode++)); do
    for((ord=1;ord<=max;ord++)); do
        if(( mode )); then
            printf "  if(func.info->flags & XC_FLAGS_HAVE_%s) {\n" ${der[ord-1]}
        fi
        for((tau=0;tau<=max;tau++)); do
            for((lapl=0;lapl<=max-tau;lapl++)); do
                for((sigma=0;sigma<=max-tau-lapl;sigma++)); do
                    for((rho=0;rho<=max-tau-lapl-sigma;rho++)); do
                        let tot=rho+sigma+lapl+tau
                        if(( tot != ord )); then
                            continue
                        fi

                        if(( ord == 1 )); then
                            name="v"
                        else
                            name=$(printf "v%i" $ord)
                        fi

                        if(( rho )); then
                            if(( rho == 1 )); then
                                name=$(printf "%srho" $name)
                            else
                                name=$(printf "%srho%i" $name $rho)
                            fi
                        fi
                        if(( sigma )); then
                            if(( sigma == 1 )); then
                                name=$(printf "%ssigma" $name)
                            else
                                name=$(printf "%ssigma%i" $name $sigma)
                            fi
                        fi
                        if(( lapl )); then
                            if(( lapl == 1 )); then
                                name=$(printf "%slapl" $name)
                            else
                                name=$(printf "%slapl%i" $name $lapl)
                            fi
                        fi
                        if(( tau )); then
                            if(( tau == 1 )); then
                                name=$(printf "%stau" $name)
                            else
                                name=$(printf "%stau%i" $name $tau)
                            fi
                        fi
                        
                        if(( mode )); then
                            printf "    %-18s = values.%s;\n" p${name} $name
                        else
                            printf "  double *%-18s = NULL;\n" p${name}
                        fi
                    done
                done
            done
        done
        if(( mode )); then
            printf "  }\n"
        fi
    done
done
