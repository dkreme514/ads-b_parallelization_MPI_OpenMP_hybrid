# High Performance Parallel Air-Traffic Analysis

C++ | MPI | OpenMP | Distributed Systems | HPC

A high-performance parallel system for analyzing large-scale ADS-B aircraft telemetry data using a hybrid MPI + OpenMP architecture. This project focuses on scalable data processing, efficient geospatial computation, and distributed sorting for real-time air traffic analysis.

Overview

This project implements a hybrid parallel pipeline that processes massive aircraft telemetry datasets to perform proximity analysis, distance computation, and global sorting across distributed nodes.

The system is designed to:

- Scale efficiently across multiple cores and nodes
- Optimize computational bottlenecks in geospatial calculations
- Reduce latency for real-time or near real-time analytics

Tech Stack
- Language: C++
- Parallelism:
    MPI (distributed memory)
    OpenMP (shared memory threading)
- Concepts:
    Hybrid parallel computing
    Distributed sorting
    Load balancing
    Strong scaling analysis

Key Features
Hybrid Parallel Pipeline
- Combined MPI for inter-node communication with OpenMP for intra-node threading
- Achieved up to 11× strong scaling speedup on 32 cores
- Efficient workload partitioning across distributed processes

Geospatial Distance Optimization

Implemented and benchmarked four distance algorithms:
- Vincenty (high accuracy, higher cost)
- Haversine (balanced accuracy/performance)
- Half-Haversine (optimized variant)
- Equirectangular approximation (fast, lower accuracy)

Result:

- Reduced proximity query latency by over 40%
- Enabled dynamic trade-offs between accuracy and performance

Distributed Sample Sort
- Built a scalable sorting pipeline using:
    MPI_Alltoallv for data redistribution
    Per-thread priority queues for local sorting
- Improved global sort efficiency by ~10% vs merge-based approaches

Performance Highlights
Metric	                  Result
Strong Scaling	          11× speedup (32 cores)
Proximity Query Latency	  ↓ 40%
Sorting Efficiency	      ↑ 10% vs baseline

System Design

The system follows a multi-stage pipeline:

1. Data Partitioning (MPI)
    - Distribute ADS-B data across nodes
2. Parallel Processing (OpenMP)
    - Multi-threaded distance computations
    - Local filtering and aggregation
3. Distributed Sorting
    - Sample-based partitioning
    - All-to-all communication
    - Local priority queue sorting
4. Aggregation & Output
    - Final global ordering and analysis

Project Structure
/src
  ├── main.cpp
  ├── mpi_pipeline.cpp
  ├── distance_algorithms.cpp
  ├── sample_sort.cpp
  ├── utils.cpp

/include
  ├── distance_algorithms.h
  ├── sample_sort.h

/data
  ├── sample_adsb_data.csv

/results
  ├── performance_logs.txt

How to Run
Prerequisites
- C++ compiler (GCC/Clang with OpenMP support)
- MPI implementation (e.g., MPICH or OpenMPI)

Build
  mpicxx -fopenmp -O3 -o air_traffic main.cpp

Execute
  mpirun -np 4 ./air_traffic

Example Use Cases
- Real-time aircraft proximity detection
- Airspace congestion analysis
- Distributed geospatial analytics
- High-performance data pipelines for streaming telemetry

Future Improvements
- GPU acceleration (CUDA/OpenACC)
- Integration with streaming frameworks (e.g., Kafka)
- Cloud deployment (AWS parallel clusters / Kubernetes)
- Real-time visualization dashboard
