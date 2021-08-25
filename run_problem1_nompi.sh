#!/bin/bash
#PJM --rsc-list "node=1"
#PJM --rsc-list "rscgrp=supercon2021"
#PJM --rsc-list "elapse=0:06:00"
#PJM -j

export OMP_NUM_THREADS=48
echo $OMP_NUM_THREADS

# do not generate empty output files
export PLE_MPI_STD_EMPTYFILE=off

date

# compile: choose FCC (for C++) or fcc (for C) and edit the source file name
#  C++
FCC -Ofast -fopenmp -Nclang -SSL2BLAMP prog.cpp
#  C
#fcc -Ofast -fopenmp -Nclang -SSL2BLAMP prog.c

# execute the program
./a.out < input1 > output1.${PJM_JOBID}


