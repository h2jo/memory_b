#!/bin/bash
gcc main.c -o main.x -lm -O4

ens0="0"
ens1="10"
mu="20"
nu="30"
T="500000"

for rho in 0.000 0.001 0.010
do
    ./main.x $ens0 $ens1 $mu $nu $rho $T
done
