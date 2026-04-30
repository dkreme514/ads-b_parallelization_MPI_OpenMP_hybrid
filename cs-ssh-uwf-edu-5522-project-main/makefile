# cs-ssh.uwf.edu --- module load mpi/openmpi-x86_64

CC        ?= gcc
MPICC     ?= mpicc
CFLAGS    ?= -std=c11 -O3 -Wall -Wextra -Wpedantic -Wformat=2 -march=native
LDFLAGS   ?= -lm
OMPFLAGS  ?= -fopenmp

# Sources
SRC_COMMON := common.c
OBJ_COMMON := common.o

SRC_SERIAL := find_closest_serial.c
BIN_SERIAL := find_closest_serial

SRC_OMP := find_closest_omp.c
BIN_OMP := find_closest_omp

SRC_MPI := find_closest_mpi.c
BIN_MPI := find_closest_mpi

BINS := $(BIN_SERIAL) $(BIN_OMP) $(BIN_MPI)

.PHONY: all clean dirs

all: $(BINS)

dirs:
	@mkdir -p input output

$(OBJ_COMMON): $(SRC_COMMON) common.h | dirs
	$(CC) $(CFLAGS) $(OMPFLAGS) -c $(SRC_COMMON) -o $(OBJ_COMMON)

$(BIN_SERIAL): $(SRC_SERIAL) $(OBJ_COMMON) | dirs
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BIN_OMP): $(SRC_OMP) $(OBJ_COMMON) | dirs
	$(CC) $(CFLAGS) $(OMPFLAGS) $^ -o $@ $(LDFLAGS)

# MPI build: use mpicc so <mpi.h> and MPI libs are wired in
find_closest_mpi.o: $(SRC_MPI) common.h | dirs
	$(MPICC) $(CFLAGS) $(OMPFLAGS) -c $(SRC_MPI) -o $@

$(BIN_MPI): find_closest_mpi.o $(OBJ_COMMON) | dirs
	$(MPICC) $(CFLAGS) $(OMPFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	@rm -f $(BINS) $(OBJ_COMMON) find_closest_mpi.o

