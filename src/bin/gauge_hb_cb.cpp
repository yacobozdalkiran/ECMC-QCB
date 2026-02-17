#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../flow/gradient_flow.h"
#include "../gauge/GaugeField.h"
#include "../heatbath/heatbath_mpi.h"
#include "../io/ildg.h"
#include "../io/io.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/HalosShift.h"
#include "../mpi/Shift.h"
#include "../observables/observables_mpi.h"

namespace fs = std::filesystem;

void generate_hb_cb(const RunParamsHbCB& rp, bool existing) {
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
    // Initalization of halos
    mpi::exchange::exchange_halos_cascade(field, geo, topo);

    // Params HB
    HbParams hp = rp.hp;

    // Shift objects
    HalosShift halo_shift(geo);

    int N_shift = rp.N_shift;
    int N_switch_eo = rp.N_switch_eo;

    // Topo
    double eps = 0.02;
    GradientFlow flow(eps, field, geo);

    // Measure vectors
    std::vector<double> plaquette;
    std::vector<double> plaquette_even;
    std::vector<double> plaquette_odd;
    plaquette_even.reserve(rp.hp.N_samples);
    plaquette_odd.reserve(rp.hp.N_samples);
    plaquette.reserve(rp.N_shift * rp.N_switch_eo * 2 * rp.hp.N_samples);

    std::vector<double> tQE_tot;
    std::vector<double> tQE_current;
    if (rp.topo) {
        tQE_tot.reserve((rp.N_shift / rp.N_shift_topo) * 3 * rp.N_rk_steps * rp.N_steps_gf);
        tQE_current.reserve(3 * rp.N_rk_steps * rp.N_steps_gf);
    }

    // Print params
    print_parameters(rp, topo);

    //==============================Heatbath Checkboard===========================

    // Thermalisation

    if (topo.rank == 0) {
        std::cout << "Thermalisation : " << rp.N_therm << " shifts\n";
    }
    for (int i = 0; i < rp.N_therm; i++) {
        for (int j = 0; j < N_switch_eo; j++) {
            // Even parity :
            if (topo.rank == 0) {
                std::cout << "(Therm) Shift : " << i << ", Switch : " << j << ", Parity : Even\n";
            }
            parity active_parity = even;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette_even = mpi::heatbathcb::samples(field, geo, topo, hp, rng, active_parity);
            // plaquette.insert(plaquette.end(), std::make_move_iterator(plaquette_even.begin()),
            // std::make_move_iterator(plaquette_even.end()));
            // Odd parity :
            if (topo.rank == 0) {
                std::cout << "(Therm) Shift : " << i << ", Switch : " << j << ", Parity : Odd\n";
            }
            active_parity = odd;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette_odd = mpi::heatbathcb::samples(field, geo, topo, hp, rng, active_parity);
            // plaquette.insert(plaquette.end(), std::make_move_iterator(plaquette_odd.begin()),
            // std::make_move_iterator(plaquette_odd.end()));
        }

        // Random shift
        mpi::shift::random_shift(field, geo, halo_shift, topo, rng);
    }

    // Sampling
    if (topo.rank == 0) {
        std::cout << "Sampling : " << rp.N_shift * rp.N_switch_eo * 2 * rp.hp.N_samples
                  << " <P> samples, " << rp.N_shift / rp.N_shift_topo << " Q samples\n";
    }

    for (int i = 0; i < N_shift; i++) {
        for (int j = 0; j < N_switch_eo; j++) {
            // Even parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Even\n";
            }
            parity active_parity = even;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette_even = mpi::heatbathcb::samples(field, geo, topo, hp, rng, active_parity);
            plaquette.insert(plaquette.end(), std::make_move_iterator(plaquette_even.begin()),
                             std::make_move_iterator(plaquette_even.end()));
            // Odd parity :
            if (topo.rank == 0) {
                std::cout << "Shift : " << i << ", Switch : " << j << ", Parity : Odd\n";
            }
            active_parity = odd;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette_odd = mpi::heatbathcb::samples(field, geo, topo, hp, rng, active_parity);
            plaquette.insert(plaquette.end(), std::make_move_iterator(plaquette_odd.begin()),
                             std::make_move_iterator(plaquette_odd.end()));
        }

        // Random shift
        mpi::shift::random_shift(field, geo, halo_shift, topo, rng);

        // Measure topo
        if (rp.topo and (i % rp.N_shift_topo == 0)) {
            tQE_current = mpi::observables::topo_charge_flowed(field, geo, flow, topo,
                                                               rp.N_steps_gf, rp.N_rk_steps);
            tQE_tot.insert(tQE_tot.end(), std::make_move_iterator(tQE_current.begin()),
                           std::make_move_iterator(tQE_current.end()));
        }
    }

    //===========================Output======================================

    if (topo.rank == 0) {
        // Write the output
        int precision = 10;
        io::save_double(plaquette, rp.run_name, precision);
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
    RunParamsHbCB params;
    bool existing = io::read_params(params, rank, argv[1]);

    // Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    generate_hb_cb(params, existing);

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
