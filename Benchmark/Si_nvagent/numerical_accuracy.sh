#!/usr/bin/env bash
# numerical_accuracy.sh — summarize the SCF iterations in a MuST main output (o_*) file.
#
# Usage:
#   ./numerical_accuracy.sh o_n0000000_Si                         # print each iteration
#   ./numerical_accuracy.sh o_n0000000_Si reference_run/o_n0000000_Si
#                                                                 # also compare to a reference
#
# For every SCF iteration it prints, in file order:
#   Fermi Energy DOS per Site, Fermi Energy (Ryd), Integrated Electron Density in VP,
#   Integrated Moment Density in VP, Interstitial rho0, Total Energy Per Atom
# An iteration is taken to end at its "Total Energy Per Atom" line.
#
# With a second (reference) file, it also prints the per-iteration difference
# |run - reference| of the total energy and the Fermi energy.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 <o_file> [reference_o_file]" >&2
    exit 1
fi

extract() {
    awk '
    BEGIN {
        n = split("Fermi Energy DOS per Site|Fermi Energy (Ryd)|Integrated Electron Density in VP|" \
                  "Integrated Moment Density in VP|Interstitial rho0|Total Energy Per Atom", keys, "|")
        it = 0; buf = ""
    }
    {
        line = $0; sub(/[ \t\r]+$/, "", line)
        for (i = 1; i <= n; i++) {
            k = keys[i]
            # line must start (after spaces) with the key, then spaces and "="
            if (index(line, k) > 0 && line ~ ("^[ \t]*" esc(k) "[ \t]*="))  {
                buf = buf line "\n"
                if (k == "Total Energy Per Atom") {
                    it++
                    printf "=================== ITERATION %02d ===================\n%s\n", it, buf
                    buf = ""
                }
                break
            }
        }
    }
    END { if (it == 0) { print "No \"Total Energy Per Atom\" lines found." > "/dev/stderr"; exit 2 } }
    function esc(s) { gsub(/[][()\\.^$*+?{}|]/, "\\\\&", s); return s }
    ' "$1"
}

# Pull "<iteration> <value>" pairs for one key (first occurrence per iteration)
series() {   # $1 = file, $2 = key
    extract "$1" | awk -v key="$2" '
        /^=+ ITERATION/ { it = $3 + 0; got = 0; next }
        !got && index($0, key) { split($0, a, "="); split(a[2], v, " "); printf "%02d %s\n", it, v[1]; got = 1 }'
}

run="$1"
[[ -f "$run" ]] || { echo "File not found: $run" >&2; exit 1; }

echo "File: $run"
echo
extract "$run"

# Final SCF energy from the matching k_ file, if present (k_ file = o_ name with k_ prefix)
kfile="$(dirname "$run")/k_$(basename "$run" | sed 's/^o_//')"
if [[ -f "$kfile" ]]; then
    echo "Final SCF energy line from $(basename "$kfile"):"
    grep -v '^[[:space:]]*$' "$kfile" | tail -n 1
    echo
fi

if [[ $# -eq 2 ]]; then
    ref="$2"
    [[ -f "$ref" ]] || { echo "Reference file not found: $ref" >&2; exit 1; }
    echo "=================== COMPARISON vs $ref ==================="
    printf "%-5s %22s %22s %12s   %18s %18s %12s\n" "Iter" "E_tot run (Ry)" "E_tot ref (Ry)" "|dE|" \
           "E_F run (Ry)" "E_F ref (Ry)" "|dE_F|"
    join <(series "$run" "Total Energy Per Atom") <(series "$ref" "Total Energy Per Atom") > /tmp/.na_e.$$
    join <(series "$run" "Fermi Energy (Ryd)")    <(series "$ref" "Fermi Energy (Ryd)")    > /tmp/.na_f.$$
    join /tmp/.na_e.$$ /tmp/.na_f.$$ | awk '
        function abs(x) { return x < 0 ? -x : x }
        { printf "%-5d %22s %22s %12.3e   %18s %18s %12.3e\n", $1, $2, $3, abs($2-$3), $4, $5, abs($4-$5) }'
    rm -f /tmp/.na_e.$$ /tmp/.na_f.$$
fi
