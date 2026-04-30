#!/usr/bin/env bash
#
# Sweep MPI ranks (-np) and OpenMP threads (OMP_NUM_THREADS)
# to study scaling of find_closest_mpi.
#
# Usage:  bash sweep_mpi_omp.sh
#   (edit vars below to match your paths / settings)

# ---- user-tunable parameters ----
BINARY="./find_closest_mpi"

DATA_FILE="input/input-data.csv"
PHASE=1               # 1=Vincenty, 2=Full Hav, 3=Half a, 4=Equirect
REF_LAT=30.4733
REF_LON=-87.1866
TOPX=25

MAX_NP=10             # max MPI ranks to test (1..MAX_NP)
MAX_THREADS=10        # max OMP threads per rank (1..MAX_THREADS)

LOG_STDOUT="scaling_runs.log"
TIMING_CSV="output/timings.csv"

# If you know the max cores per node, you can set this to avoid oversubscription.
# Set to 0 to disable check.
TOTAL_CORES_PER_NODE=0   # e.g. 32; 0 = don't enforce

# ---------------------------------

echo "Sweeping np=1..$MAX_NP, OMP_NUM_THREADS=1..$MAX_THREADS"
echo "Binary:     $BINARY"
echo "Data file:  $DATA_FILE"
echo "Phase:      $PHASE"
echo "Ref:        ($REF_LAT, $REF_LON), TopX=$TOPX"
echo

# Make sure output dir exists
mkdir -p output

# Clear old logs
rm -f "$LOG_STDOUT"
rm -f "$TIMING_CSV"

echo "# Scaling sweep $(date)"        | tee -a "$LOG_STDOUT"
echo "# BINARY=$BINARY"              | tee -a "$LOG_STDOUT"
echo "# DATA_FILE=$DATA_FILE"        | tee -a "$LOG_STDOUT"
echo "# PHASE=$PHASE"                | tee -a "$LOG_STDOUT"
echo "# REF_LAT=$REF_LAT"            | tee -a "$LOG_STDOUT"
echo "# REF_LON=$REF_LON"            | tee -a "$LOG_STDOUT"
echo "# TOPX=$TOPX"                  | tee -a "$LOG_STDOUT"
echo "# MAX_NP=$MAX_NP"              | tee -a "$LOG_STDOUT"
echo "# MAX_THREADS=$MAX_THREADS"    | tee -a "$LOG_STDOUT"
echo >> "$LOG_STDOUT"

# Ensure OpenMP doesn’t change thread counts dynamically
export OMP_DYNAMIC=FALSE

for np in $(seq 1 "$MAX_NP"); do
  for threads in $(seq 1 "$MAX_THREADS"); do

    # Optional: skip combinations that exceed a known core limit
    if [ "$TOTAL_CORES_PER_NODE" -gt 0 ]; then
      total=$((np * threads))
      if [ "$total" -gt "$TOTAL_CORES_PER_NODE" ]; then
        echo "Skipping np=$np, threads=$threads (total=$total > $TOTAL_CORES_PER_NODE cores)" \
          | tee -a "$LOG_STDOUT"
        continue
      fi
    fi

    echo "========================================================" | tee -a "$LOG_STDOUT"
    echo "Run: np=$np, OMP_NUM_THREADS=$threads"                    | tee -a "$LOG_STDOUT"
    echo "--------------------------------------------------------" | tee -a "$LOG_STDOUT"

    export OMP_NUM_THREADS="$threads"

    # You can swap mpirun/mpiexec to match your environment
    mpirun -np "$np" \
      --bind-to none \
      "$BINARY" \
      --phase="$PHASE" \
      "$REF_LAT" "$REF_LON" "$TOPX" \
      --data-file "$DATA_FILE" \
      2>&1 | tee -a "$LOG_STDOUT"

    echo                                                            | tee -a "$LOG_STDOUT"
  done
done

echo "Done. Stdout in $LOG_STDOUT, timings aggregated in $TIMING_CSV"

