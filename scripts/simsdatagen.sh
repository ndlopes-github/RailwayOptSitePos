#!/bin/bash
# This script is used to generate the simulated instances data

echo "loop start"
for c in {1..3} #{1..11} 
do
    for i in {1..3} #{1..10}
    do
        command time -v  julia ./sim_instances_generator.jl $((2**c)) $((i)) &
    done
    wait
done
echo "loop end"
