#!/bin/bash


mpifort -c threeD_FFTW_CHEB.f90 -I/usr/local/include -L/usr/local/lib -lfftw3_mpi -lfftw3 -llapack -lblas -O3 -march=native
mpifort -c QG_plane.f90 -I/usr/local/include -L/usr/local/lib -llapack -lblas -O3
mpifort -o prova QG_plane.o threeD_FFTW_CHEB.o -I/usr/local/include -L/usr/local/lib -lfftw3_mpi -lfftw3 -llapack -lblas -O3 -march=native
mpifort -c QG_plane_front.f90 -I/usr/local/include -L/usr/local/lib -llapack -lblas -O3
mpifort -o prova_front QG_plane_front.o threeD_FFTW_CHEB.o -I/usr/local/include -L/usr/local/lib -lfftw3_mpi -lfftw3 -llapack -lblas -O3 -march=native

echo "Compilation completed"
