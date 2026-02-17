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
#include "../mpi/MpiTopology.h"
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
        read_ildg_clime("data/"+rp.run_name, field, geo, topo);
        fs::path state_path =
            fs::path("data/"+rp.run_name+"_seed") / (rp.run_name + "_seed" + std::to_string(topo.rank) + ".txt");
        std::ifstream ifs(state_path);
        if (ifs.is_open()) {
            ifs >> rng;  // On ignore la seed initiale, on reprend l'état exactement
        } else {
            std::cerr << "Could not open " << state_path << "\n";
        }
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
        io::save_topo(tQE_tot, rp.run_name, precision);
        io::save_params(rp, rp.run_name);
    }
    //Save seeds
    io::save_seed(rng, rp.run_name, topo);
    // Save conf
    save_ildg_clime("data/" + rp.run_name, field, geo, topo);
}

// Reads the parameters of input file into RunParams struct
// If run_name.ildg is found returns true
bool read_params(RunParamsHbCB& params, int rank, const std::string& input) {
    if (rank == 0) {
        try {
            io::load_params(input, params);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input : " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // 1. Diffusion des paramètres numériques (Lattice + Run + Topo)
    // On diffuse les blocs un par un pour plus de clarté et de sécurité
    MPI_Bcast(&params.L_core, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.n_core_dims, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.cold_start, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.N_switch_eo, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.N_shift, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.N_therm, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.topo, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.N_shift_topo, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.N_steps_gf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.N_rk_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&params.hp.beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.hp.N_samples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.hp.N_hits, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.hp.N_sweeps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&params.hp.N_therm, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 3. Diffusion de la std::string (run_name)
    int name_len;
    if (rank == 0) {
        name_len = static_cast<int>(params.run_name.size());
    }

    // On envoie d'abord la taille de la string
    MPI_Bcast(&name_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Les esclaves préparent leur mémoire
    if (rank != 0) {
        params.run_name.resize(name_len);
    }

    // On envoie le contenu de la string (le buffer interne)
    if (name_len > 0) {
        MPI_Bcast(&params.run_name[0], name_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    bool local_existing =
        fs::exists("data/" + params.run_name) and
        fs::exists("data/" + params.run_name + "_plaquette.txt") and
        fs::exists("data/" + params.run_name + "_topo.txt") and
        fs::exists("data/" + params.run_name + "_seed/" + params.run_name + "_seed" + std::to_string(rank) + ".txt");
    bool global_existing;
    // On vérifie que TOUS les processus ont trouvé leurs fichiers
    MPI_Allreduce(&local_existing, &global_existing, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    if (global_existing) {
        if (rank == 0) std::cout << "All ranks found existing configuration. Resuming...\n";
        return true;
    } else {
        return false;
    }
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
    bool existing = read_params(params, rank, argv[1]);

    // Measuring time
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    generate_hb_cb(params, existing);

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
