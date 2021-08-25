#include <omp.h>
#include "mpi.h"
#include "cssl.h"
#include "sc21.h"
#include <fstream>
#include <iostream>
using namespace std;

int main(int argc, char **argv) {
    int myid, n_procs;
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&n_procs); 
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    cout << myid << " " << n_procs << endl;

    MPI_Finalize();
    return 0;
}
