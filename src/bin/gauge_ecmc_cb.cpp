#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "Projet ECMC_QCB initialisé !" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
