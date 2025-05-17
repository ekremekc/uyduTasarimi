# python3 mesh.py -nopopup -thickness 0.001 && mpirun -np 8 python3 main.py -thickness 0.001


#!/bin/bash

start=0.0001
end=0.001
step=0.0001

# Generate the array using seq and printf
thickness_values=($(seq $start $step $end | xargs printf "%.4f\n"))

for t in "${thickness_values[@]}"; do
    echo "Running for thickness = $t"
    python3 mesh.py -nopopup -thickness "$t" && mpirun -np 8 python3 main.py -thickness "$t"
done
