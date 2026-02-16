#include <iostream>
#include <mpi.h>

#include "../gauge/GaugeField.h"
#include "../mpi/HalosExchange.h"
#include "../observables/observables_mpi.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

    int L = 4;
    GeometryCB geo(L);
    GaugeField field(geo);

    int n_core_dims = 2;
    mpi::MpiTopology topo(n_core_dims);

    std::random_device rd;
    std::mt19937_64 rng(rd()+topo.rank);

    field.hot_start(geo, rng);

    double p = mpi::observables::mean_plaquette_global(field, geo, topo);
    if (topo.rank == 0){
        std::cout << "P = " << p << "\n";
    }

    MPI_Finalize();
    return 0;
}
