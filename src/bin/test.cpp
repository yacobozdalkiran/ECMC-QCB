
#include <iomanip>
#include <iostream>

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

    // MPI
    int n_core_dims = 2;
    mpi::MpiTopology topo(n_core_dims);

    // Lattice
    int L_core = 8;
    GeometryCB geo(L_core);
    GaugeField field(geo);

    // Flow
    GradientFlow flow(0.02, field, geo);

    // Load
    std::string filename = "P3r9a";
    std::string filedir = "data";
    read_ildg_clime(filename, filedir, field, geo, topo);

    // Checks
    double p = mpi::observables::mean_plaquette_global(field, geo, topo);
    if (topo.rank == 0) {
        std::cout << "P = " << p << "\n";
        std::cout << "U_x(0) " << "\n";
        std::cout << field.view_link_const(geo.index(1, 1, 1, 1), 0) << "\n";
        std::cout << "U_t(0) " << "\n";
        std::cout << field.view_link_const(geo.index(1, 1, 1, 1), 3) << "\n";
    }
    // Unitarity
    double error_max = 0.0;
    double error_mean = 0.0;
    double error_mean_sq = 0.0;
    double count = 0;

    field.project_field_su3(geo);
    for (int t = 1; t < geo.L_ext - 1; t++) {
        for (int z = 1; z < geo.L_ext - 1; z++) {
            for (int y = 1; y < geo.L_ext - 1; y++) {
                for (int x = 1; x < geo.L_ext - 1; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        SU3 U = field.view_link_const(site, mu);
                        double error_link =
                            std::sqrt(((U * U.adjoint() - SU3::Identity()).adjoint() *
                                       (U * U.adjoint() - SU3::Identity()))
                                          .trace()
                                          .real());
                        if (error_max < error_link) error_max = error_link;
                        error_mean += error_link;
                        error_mean_sq += error_link * error_link;
                        count += 1;
                    }
                }
            }
        }
    }
    error_mean /= count;
    error_mean_sq /= count;

    double error_var = error_mean_sq - error_mean * error_mean;
    double error_max_global = 0.0;
    double error_mean_global = 0.0;
    double var_global = 0.0;
    double error_std_global = 0.0;

    MPI_Reduce(&error_mean, &error_mean_global, 1, MPI_DOUBLE, MPI_SUM, 0, topo.cart_comm);
    MPI_Reduce(&error_var, &var_global, 1, MPI_DOUBLE, MPI_SUM, 0, topo.cart_comm);
    MPI_Reduce(&error_max, &error_max_global, 1, MPI_DOUBLE, MPI_MAX, 0, topo.cart_comm);

    if (topo.rank == 0) {
        error_std_global = std::sqrt(var_global)/(double)topo.size;
        error_mean_global /= (double)topo.size;
        std::cout << "Error max : " << error_max_global << "\n";
        std::cout << "Error mean : " << error_mean_global << "\n";
        std::cout << "Error std : " << error_std_global << "\n";
    }

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
