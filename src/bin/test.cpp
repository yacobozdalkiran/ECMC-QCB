
#include <iostream>
#include <iomanip>

#include "../flow/gradient_flow.h"
#include "../gauge/GaugeField.h"
#include "../io/ildg.h"
#include "../io/io.h"
#include "../mpi/MpiTopology.h"
#include "../observables/observables_mpi.h"


int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    //MPI
    int n_core_dims = 2;
    mpi::MpiTopology topo(n_core_dims);

    //Lattice
    int L_core = 4;
    GeometryCB geo(L_core);
    GaugeField field(geo);

    //Flow
    GradientFlow flow(0.02, field, geo);

    //Load
    std::string filename="data/conf.ildg";
    read_ildg_clime(filename, field, geo, topo);
    
    //Checks
    double p = mpi::observables::mean_plaquette_global(field, geo, topo);
    if (topo.rank == 0){
        std::cout << "P = " << p << "\n";
        std::cout << "U_x(0) " << "\n";
        std::cout << field.view_link_const(geo.index(1,1,1,1),0) << "\n";
        std::cout << "U_t(0) " << "\n";
        std::cout << field.view_link_const(geo.index(1,1,1,1),3) << "\n";
    }
    mpi::observables::topo_charge_flowed(field, geo, flow, topo, 11, 40);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (topo.rank == 0) {
        double total_time = end_time - start_time;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\n==========================================" << std::endl;
        std::cout << " Total execution time : " << total_time << " seconds" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
    // End MPI
    MPI_Finalize();
}
