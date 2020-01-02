#!/bin/bash

# g++-8 -O2 -I ../../../include renormalizedAtom.cpp -o renormalizedAtom
# g++-8 -O2 -I ../../../include renormalizedAtomNonRelativistic.cpp -o renormalizedAtomNonRelativistic
g++-8 -O2 -I ../../../include finiteDifferenceAtom.cpp -o finiteDifferenceAtom -framework Accelerate
g++-8 -O2 -I ../../../include finiteDifferenceAtomSpinPolarized.cpp -o finiteDifferenceAtomSpinPolarized -framework Accelerate

