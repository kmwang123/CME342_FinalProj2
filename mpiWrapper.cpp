#include <mpi.h>
#include "mpiWrapper.hpp"

using namespace std;

MPI_Wrapper::MPI_Wrapper(void) {
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  comm = MPI_COMM_WORLD;
}
