#include <mpi.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../ecmc/ecmc_mpi_cb.h"
#include "../flow/gradient_flow.h"
#include "../gauge/GaugeField.h"
#include "../io/ildg.h"
#include "../io/io.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/HalosShift.h"
#include "../mpi/Shift.h"
#include "../observables/observables_mpi.h"

namespace fs = std::filesystem;

void generate_ecmc_cb(const RunParamsECB& rp, bool existing) {
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

    if (existing) {
        read_ildg_clime("data/" + rp.run_name + "/" + rp.run_name, field, geo, topo);
        fs::path state_path = fs::path("data/" + rp.run_name + "/" + rp.run_name + "_seed") /
                              (rp.run_name + "_seed" + std::to_string(topo.rank) + ".txt");
        std::ifstream ifs(state_path);
        if (ifs.is_open()) {
            ifs >> rng;  // On ignore la seed initiale, on reprend l'état exactement
        } else {
            std::cerr << "Could not open " << state_path << "\n";
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Initalization of halos for ECMC
    mpi::exchange::exchange_halos_cascade(field, geo, topo);

    // Params ECMC
    ECMCParams ep = rp.ecmc_params;

    // Shift objects
    HalosShift halo_shift(geo);

    int N_shift = rp.N_shift;
    int N_switch_eo = rp.N_switch_eo;

    // Topo
    double eps = 0.02;
    GradientFlow flow(eps, field, geo);

    // Measure vectors
    std::vector<double> plaquette;
    plaquette.reserve(rp.N_shift);

    std::vector<double> tQE_tot;
    std::vector<double> tQE_current;
    if (rp.topo) {
        tQE_tot.reserve((rp.N_shift / rp.N_shift_topo) * 3 * rp.N_rk_steps * rp.N_steps_gf);
        tQE_current.reserve(3 * rp.N_rk_steps * rp.N_steps_gf);
    }

    // Print params
    print_parameters(rp, topo);

    //==============================ECMC Checkboard===========================
    // Thermalisation

    if (topo.rank == 0) {
        std::cout << "Thermalisation : " << rp.N_therm << " shifts\n";
    }

    for (int i = 0; i < rp.N_therm; i++) {
        if (topo.rank == 0){
            std::cout << "==========" << "(Therm) Sample " << i << "==========\n";
        }
        for (int j = 0; j < N_switch_eo; j++) {
            // Even parity :
            if (topo.rank == 0) {
                std::cout << "(Therm) Shift : " << i << ", Switch : " << j << ", Parity : Even\n";
            }
            parity active_parity = even;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            mpi::ecmccb::sample(field, geo, ep, rng, topo, active_parity);

            // Odd parity :
            if (topo.rank == 0) {
                std::cout << "(Therm) Shift : " << i << ", Switch : " << j << ", Parity : Odd\n";
            }
            active_parity = odd;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            mpi::ecmccb::sample(field, geo, ep, rng, topo, active_parity);
        }

        // Plaquette measure (not saved for thermalization)
        double p = mpi::observables::mean_plaquette_global(field, geo, topo);
        if (topo.rank == 0) {
            std::cout << "(Therm) Sample " << i << ", <P> = " << p << " ";
            std::cout << "\n";
        }
        // Random shift
        mpi::shift::random_shift(field, geo, halo_shift, topo, rng);
    }

    // Sampling
    if (topo.rank == 0) {
        std::cout << "Sampling : " << rp.N_shift << " <P> samples, " << rp.N_shift / rp.N_shift_topo
                  << " Q samples\n";
    }

    for (int i = 0; i < N_shift; i++) {
        if (topo.rank == 0){
            std::cout << "=============" << "Sample " << i << "=============\n";
        }
        for (int j = 0; j < N_switch_eo; j++) {
            // Even parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Even\n";
            }
            parity active_parity = even;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            mpi::ecmccb::sample(field, geo, ep, rng, topo, active_parity);

            // Odd parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Odd\n";
            }
            active_parity = odd;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            mpi::ecmccb::sample(field, geo, ep, rng, topo, active_parity);
        }

        // Plaquette measure
        double p = mpi::observables::mean_plaquette_global(field, geo, topo);
        if (topo.rank == 0) {
            std::cout << "Sample " << i << ", <P> = " << p << " ";
            std::cout << "\n";
        }
        plaquette.emplace_back(p);
        // Measure topo
        if (rp.topo and (i % rp.N_shift_topo == 0)) {
            if (topo.rank == 0){
                std::cout << "Measuring Q sample " << i/rp.N_shift_topo << "\n";
            }
            tQE_current = mpi::observables::topo_charge_flowed(field, geo, flow, topo,
                                                               rp.N_steps_gf, rp.N_rk_steps);
            tQE_tot.insert(tQE_tot.end(), std::make_move_iterator(tQE_current.begin()),
                           std::make_move_iterator(tQE_current.end()));
        }
        
        // Random shift
        mpi::shift::random_shift(field, geo, halo_shift, topo, rng);
    }

    //===========================Output======================================

    if (topo.rank == 0) {
        // Write the output
        int precision = 10;
        io::save_plaquette(plaquette, rp.run_name, precision);
        if (rp.topo) {
            io::save_topo(tQE_tot, rp.run_name, precision);
        }
        io::save_params(rp, rp.run_name);
    }
    // Save seeds
    io::save_seed(rng, rp.run_name, topo);
    // Save conf
    save_ildg_clime("data/" + rp.run_name + "/" + rp.run_name, field, geo, topo);
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
    bool existing = io::read_params(params, rank, argv[1]);

    // Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    generate_ecmc_cb(params, existing);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double total_time = end_time - start_time;
        std::cout << std::fixed << std::setprecision(4);
        print_time(total_time);
    }
    // End MPI
    MPI_Finalize();
}
