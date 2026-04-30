#!/bin/bash
#SBATCH -N 1                  # single node
#SBATCH -n 1                  # single MPI rank (size=1)
#SBATCH --cpus-per-task=128   # allow up to 128 OpenMP threads
#SBATCH -t 00:15:00
#SBATCH -p RM
#SBATCH -J weak_omp
#SBATCH -o output/batch_weak_omp.%j.out
#SBATCH -e output/batch_weak_omp.%j.err

module load gcc
make -B
mkdir -p output
lscpu

# Weak scaling parameters for OpenMP-only
#ROWS_PER_THREAD=5000000   # each thread keeps this many rows
#echo "Weak scaling of OpenMP on single node (each thread keeps ${ROWS_PER_THREAD} rows, cols=${COLS})" 

for T in   1 2 4 8 12 16 20 24 28 32 36 40 44 48 56 64 72  ; do
    export OMP_NUM_THREADS=${T}
#    N=$((ROWS_PER_THREAD * T))
    echo
    echo ">>> Running ${T} threads " 
    ./find_closest_omp --phase=2 --data-file 'input/states_2022-06-27*.csv.gz'
done


#./find_closest_omp \
#  --phase=4 \
#  --data-file 'input/states_2022-*.csv.gz' \
#  --use-rows 3000000


echo "ALL COMPLETE"

