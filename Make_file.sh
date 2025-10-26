#!/bin/bash
# Script to run create_profiles_and_Ts.py and compile Fortran programs

set -e  # stop on first error

echo "Running Python pre-processing script..."
python3 create_profiles_and_Ts.py
echo "Python script completed successfully."

echo "Starting Fortran compilation..."

# Compile threeD_FFTW_CHEB
mpifort -c threeD_FFTW_CHEB.f90 -I/usr/local/include -L/usr/local/lib \
-lfftw3_mpi -lfftw3 -llapack -lblas -O3 -march=native

# Compile QG_plane_front
mpifort -c QG_plane_front.f90 -I/usr/local/include -L/usr/local/lib \
-llapack -lblas -O3

# Link everything into executable "prova"
mpifort -o QG_planar_front_executable QG_plane_front.o threeD_FFTW_CHEB.o -I/usr/local/include -L/usr/local/lib \
-lfftw3_mpi -lfftw3 -llapack -lblas -O3 -march=native

echo "Compilation completed successfully."
