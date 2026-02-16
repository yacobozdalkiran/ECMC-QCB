#include <mpi.h>

#include <iostream>

#include "../ecmc/ecmc_mpi_cb.h"
#include "../flow/gradient_flow.h"
#include "../gauge/GaugeField.h"
#include "../io/io.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/HalosShift.h"
#include "../mpi/Shift.h"
#include "../observables/observables_mpi.h"

void print_parameters(const RunParamsECB& rp, const mpi::MpiTopology& topo) {
    if (topo.rank == 0) {
        std::cout << "==========================================" << std::endl;
        std::cout << "ECMC - Checkboard" << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "Total lattice size : " << rp.L_core * rp.n_core_dims << "^4\n";
        std::cout << "Local lattice size : " << rp.L_core << "^4\n";
        std::cout << "Beta : " << rp.ecmc_params.beta << "\n";
        std::cout << "Total number of shifts : " << rp.N_shift << "\n";
        std::cout << "Number of e/o switchs per shift : " << rp.N_switch_eo << "\n";
        std::cout << "Number of samples per checkboard step : " << rp.ecmc_params.N_samples << "\n";
        std::cout << "Total number of samples : "
                  << 2 * rp.N_switch_eo * rp.ecmc_params.N_samples * rp.N_shift << "\n";
        std::cout << "Seed : " << rp.seed << "\n";
        std::cout << "==========================================" << std::endl;
    }
}

void generate_ecmc_cb(const RunParamsECB& rp) {
    //========================Objects initialization====================
    // MPI
    int n_core_dims = rp.n_core_dims;
    mpi::MpiTopology topo(n_core_dims);

    // Lattice creation + RNG
    int L = rp.L_core;
    GeometryCB geo(L);
    GaugeField field(geo);
    std::mt19937_64 rng(rp.seed + topo.rank);
    if (!rp.cold_start) {
        field.hot_start(geo, rng);
    }

    // Initalization of halos for ECMC
    mpi::exchange::exchange_halos_cascade(field, geo, topo);

    // Params ECMC
    ECMCParams ep = rp.ecmc_params;

    // Shift objects
    HalosShift halo_shift(geo);

    int N_shift = rp.N_shift;
    int N_switch_eo = rp.N_switch_eo;

    // Print params
    print_parameters(rp, topo);

    // Measures
    std::vector<std::vector<std::vector<std::vector<double>>>> plaquette(
        rp.N_shift,
        std::vector<std::vector<std::vector<double>>>(
            rp.N_switch_eo, std::vector<std::vector<double>>(
                                2, std::vector<double>(rp.ecmc_params.N_samples, 0.0))));

    //==============================ECMC Checkboard===========================
    for (int i = 0; i < N_shift; i++) {
        for (int j = 0; j < N_switch_eo; j++) {
            // Even parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Even\n";
            }
            parity active_parity = even;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette[i][j][0] =
                mpi::ecmccb::samples_improved(field, geo, ep, rng, topo, active_parity);

            // Odd parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Odd\n";
            }
            active_parity = odd;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette[i][j][1] =
                mpi::ecmccb::samples_improved(field, geo, ep, rng, topo, active_parity);
        }
        // Random shift
        mpi::shift::random_shift(field, geo, halo_shift, topo, rng);
    }

    //===========================Gradient flow test==========================

    double eps = 0.02;
    GradientFlow flow(eps, field, geo);
    int precision = 6;
    int precision_t = 3;
    for (double time = 0.0; time < 5.0; time += eps) {
        flow.rk3_step(topo);
        auto qe = mpi::observables::topo_q_e_clover_global(flow.field_c, geo, topo);
        if (topo.rank == 0)
            std::cout << "t = " << io::format_double(time, precision_t)
                      << ", Q = " << io::format_double(qe.first, precision)
                      << ", t²E = " << io::format_double(time * time * qe.second, precision)
                      << "\n";
    }

    //===========================Output======================================

    // Flatten the vector
    if (topo.rank == 0) {
        // Flatten the plaquette vector
        std::vector<double> plaquette_flat(rp.N_shift * rp.N_switch_eo * 2 *
                                           rp.ecmc_params.N_samples);
        for (int i = 0; i < rp.N_shift; i++) {
            for (int j = 0; j < rp.N_switch_eo; j++) {
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < rp.ecmc_params.N_samples; l++) {
                        plaquette_flat[((i * rp.N_switch_eo + j) * 2 + k) *
                                           rp.ecmc_params.N_samples +
                                       l] = plaquette[i][j][k][l];
                    }
                }
            }
        }
        // Write the output
        int precision_filename = 1;
        std::string filename = "EMQCB_" + std::to_string(L * n_core_dims) + "b" +
                               io::format_double(ep.beta, precision_filename) + "Ns" +
                               std::to_string(rp.N_shift) + "Nsw" + std::to_string(rp.N_switch_eo) +
                               "Np" + std::to_string(ep.N_samples) + "c" +
                               std::to_string(rp.cold_start) + "ts" +
                               io::format_double(ep.param_theta_sample, precision_filename) + "tr" +
                               io::format_double(ep.param_theta_refresh, precision_filename);
        int precision = 10;
        io::save_double(plaquette_flat, filename, precision);
    }
}

// Reads the parameters of input file into RunParams struct
void read_params(RunParamsECB& params, int rank, const std::string& input) {
    if (rank == 0) {
        try {
            io::load_params(input, params);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input : " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // Synchronizing input parameters accross all nodes
    MPI_Bcast(&params, sizeof(RunParamsECB), MPI_BYTE, 0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Trying to read input
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <input_file.txt>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Charging the parameters of the run
    RunParamsECB params;
    read_params(params, rank, argv[1]);

    // Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    generate_ecmc_cb(params);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double total_time = end_time - start_time;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\n==========================================" << std::endl;
        std::cout << " Total execution time : " << total_time << " seconds" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
    // End MPI
    MPI_Finalize();
}
