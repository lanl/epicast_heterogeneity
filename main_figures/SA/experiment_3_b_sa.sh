#!/bin/bash 

for i in $(seq 1 10);
do  
    j=$((i + 112243))
    echo "$j" 
    mpirun -n 64 ./epicast -p "rng_seed = $j" ../../model_3_fl_sa.toml
done


