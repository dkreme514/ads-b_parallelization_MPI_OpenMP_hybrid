#!/bin/bash
#SBATCH -N 1                    # 1 node
#SBATCH --ntasks-per-node=5     # up to 5 MPI ranks for the loop
#SBATCH --cpus-per-task=24      # 24 OpenMP threads per rank
#SBATCH -t 00:10:00             # walltime hh:mm:ss
#SBATCH -p RM                   # partition/queue
#SBATCH -J weak_mpi_job
#SBATCH -o out/weak_mpi.%j.out
#SBATCH -e out/weak_mpi.%j.err

# Basic info about the node
lscpu | grep -e "CPU(s)" -e "Thread" -e "Core" -e "Socket" -e "NUMA" -e "Model"

module load gcc
module load openmpi

make -B

mkdir -p out

# Hybrid: 1 MPI rank gets 24 OpenMP threads
export OMP_NUM_THREADS=16

# Weak-scaling setup:
#   - Each rank always handles 30000 rows x 30000 columns
#   - Global rows N = 30000 * P
BASE_ROWS_PER_RANK=60000
BASE_COLS=60000

echo "Weak scaling of MPI + OpenMP" 
echo "Each rank: ${BASE_ROWS_PER_RANK} x ${BASE_COLS} local work" 
echo 

# P = number of MPI processes (you can adjust this list)
for P in 1 2 3 4 5; do
    N=$((BASE_ROWS_PER_RANK * P))
    M=$BASE_COLS
    echo "$N $M" > input.txt
    echo " $P processors on ${N} x ${M}, rows per-rank = ${BASE_ROWS_PER_RANK} " 
    mpirun --bind-to none -np "$P"  ./hw3 "$P" 
    echo "" 
done

echo "ALL COMPLETE" 

