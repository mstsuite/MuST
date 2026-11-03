#!/bin/bash

# 1. Check if an argument was provided
if [[ $# -ne 1 ]]; then
    echo "Error: Please provide exactly one o_* file name." >&2
    echo "Usage: $0 <filename>" >&2
    exit 1
fi

# 2. Check if the file actually exists and is readable
if [[ ! -f "$1" || ! -r "$1" ]]; then
    echo "Error: File '$1' does not exist or is not readable." >&2
    exit 1
fi

# 3. Stream the file in order, printing matching lines instantly
#awk '/Fermi Energy DOS per Site/ || 
#     /Integrated Electron Density in VP/ || 
#     /Integrated Moment Density in VP/ || 
#     /Fermi Energy \(Ryd\)/ || 
#     /Interstitial rho0/ || 
#     /Total Energy Per Atom/' "$1"

awk '
BEGIN { iteration = 0 }

/Fermi Energy DOS per Site|Integrated Electron Density in VP|Integrated Moment Density in VP|Fermi Energy \(Ryd\)|Interstitial rho0|Total Energy Per Atom/ {

    # Clean prefixes
    sub(/.*::[ \t]*/, "")

    # Split into Key and Value
    split($0, parts, "=")
    key = parts[1]; val = parts[2]
    gsub(/^[ \t]+|[ \t]+$/, "", key)
    gsub(/^[ \t]+|[ \t]+$/, "", val)

    # If we hit the first property, increment and print the header
    if (key == "Fermi Energy DOS per Site") {
        iteration++
        printf "\n=================== ITERATION %02d ===================\n", iteration
    }

    # Print the aligned data
    printf "  %-35s =   %s\n", key, val
}' "$1"


#grep -E "Fermi Energy DOS per Site|Integrated Electron Density in VP|Integrated Moment Density in VP|Fermi Energy \(Ryd\)|Interstitial rho0|Total Energy Per Atom" "$1"

#grep "Fermi Energy DOS per Site" $1
#grep "Integrated Electron Density in VP" $1
#grep "Integrated Moment Density in VP" $1
#grep "Fermi Energy (Ryd)" $1
#grep "Interstitial rho0" $1
#grep "Total Energy Per Atom" $1
