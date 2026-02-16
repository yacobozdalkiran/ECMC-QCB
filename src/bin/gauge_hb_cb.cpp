#include <iostream>

#include "../flow/gradient_flow.h"
#include "../gauge/GaugeField.h"
#include "../heatbath/heatbath_mpi.h"
#include "../io/io.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/HalosShift.h"
#include "../mpi/MpiTopology.h"
#include "../mpi/Shift.h"
#include "../observables/observables_mpi.h"

void print_parameters(const RunParamsHbCB& rp, const mpi::MpiTopology& topo) {
    if (topo.rank == 0) {
        std::cout << "==========================================" << std::endl;
        std::cout << "Heatbath - Checkboard" << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "Total lattice size : " << rp.L_core * rp.n_core_dims << "^4\n";
        std::cout << "Local lattice size : " << rp.L_core << "^4\n";
        std::cout << "Beta : " << rp.hp.beta << "\n";
        std::cout << "Total number of shifts : " << rp.N_shift << "\n";
        std::cout << "Number of e/o switchs per shift : " << rp.N_switch_eo << "\n";
        std::cout << "Number of sweeps : " << rp.hp.N_sweeps << "\n";
        std::cout << "Number of hits : " << rp.hp.N_hits << "\n";
        std::cout << "Number of samples per checkboard step : " << rp.hp.N_samples << "\n";
        std::cout << "Total number of samples : "
                  << 2 * rp.N_switch_eo * rp.hp.N_samples * rp.N_shift << "\n";
        std::cout << "Seed : " << rp.seed << "\n";
        std::cout << "==========================================" << std::endl;
    }
}

void generate_hb_cb(const RunParamsHbCB& rp) {
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

    // Initalization of halos
    mpi::exchange::exchange_halos_cascade(field, geo, topo);

    // Params HB
    HbParams hp = rp.hp;

    // Shift objects
    HalosShift halo_shift(geo);

    int N_shift = rp.N_shift;
    int N_switch_eo = rp.N_switch_eo;

    // Print params
    print_parameters(rp, topo);

    // Measures
    std::vector<std::vector<std::vector<std::vector<double>>>> plaquette(
        rp.N_shift, std::vector<std::vector<std::vector<double>>>(
                        rp.N_switch_eo, std::vector<std::vector<double>>(
                                            2, std::vector<double>(rp.hp.N_samples, 0.0))));

    //==============================Heatbath Checkboard===========================

    for (int i = 0; i < N_shift; i++) {
        for (int j = 0; j < N_switch_eo; j++) {
            // Even parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Even\n";
            }
            parity active_parity = even;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette[i][j][0] = mpi::heatbathcb::samples(field, geo, topo, hp, rng, active_parity);

            // Odd parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Odd\n";
            }
            active_parity = odd;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette[i][j][1] = mpi::heatbathcb::samples(field, geo, topo, hp, rng, active_parity);
        }

        // Random shift
        mpi::shift::random_shift(field, geo, halo_shift, topo, rng);
    }

    //===========================Gradient flow test==========================

    double eps = 0.02;
    GradientFlow flow(eps, field, geo);
    mpi::observables::topo_charge_flowed(field, geo, flow, topo);
    //===========================Output======================================

    // Flatten the vector
    if (topo.rank == 0) {
        // Flatten the plaquette vector
        std::vector<double> plaquette_flat(rp.N_shift * rp.N_switch_eo * 2 * rp.hp.N_samples);
        for (int i = 0; i < rp.N_shift; i++) {
            for (int j = 0; j < rp.N_switch_eo; j++) {
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < rp.hp.N_samples; l++) {
                        plaquette_flat[((i * rp.N_switch_eo + j) * 2 + k) * rp.hp.N_samples + l] =
                            plaquette[i][j][k][l];
                    }
                }
            }
        }
        // Write the output
        int precision_filename = 1;
        std::string filename = "HBQCB_" + std::to_string(L * n_core_dims) + "b" +
                               io::format_double(hp.beta, precision_filename) + "Ns" +
                               std::to_string(rp.N_shift) + "Nsw" + std::to_string(rp.N_switch_eo) +
                               "Np" + std::to_string(hp.N_samples) + "c" +
                               std::to_string(rp.cold_start) + "Nswp" +
                               std::to_string(hp.N_sweeps) + "Nh" + std::to_string(hp.N_hits);
        int precision = 10;
        io::save_double(plaquette_flat, filename, precision);
    }
}

// Reads the parameters of input file into RunParams struct
void read_params(RunParamsHbCB& params, int rank, const std::string& input) {
    if (rank == 0) {
        try {
            io::load_params(input, params);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input : " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // Synchronizing input parameters accross all nodes
    MPI_Bcast(&params, sizeof(RunParamsHbCB), MPI_BYTE, 0, MPI_COMM_WORLD);
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
    RunParamsHbCB params;
    read_params(params, rank, argv[1]);

    // Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    generate_hb_cb(params);

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
