mpiFCC -Ofast -fopenmp -Nclang -SSL2BLAMP solve.cpp
mpiexec --stdin input.txt -np 12 ./a.out
