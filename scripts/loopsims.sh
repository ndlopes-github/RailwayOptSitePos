#!/bin/bash

# This script runs all the instances used
# in the Simulated Data section.
# The required to these instances must be
# generated with the simsdatagen.sh.

echo "loop start"
for c in {1..11}  
do
    for i in {1..10}
    do
       command time -v  julia ./model_sim.jl $c $i 
    done
done
echo "loop end"
