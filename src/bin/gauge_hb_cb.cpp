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

    // Topo
    double eps = 0.02;
    GradientFlow flow(eps, field, geo);

    // Measures

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

    //Thermalisation

    if (topo.rank == 0){
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
            //plaquette.insert(plaquette.end(), std::make_move_iterator(plaquette_even.begin()),
                             std::make_move_iterator(plaquette_even.end()));
            // Odd parity :
            if (topo.rank == 0) {
                std::cout << "(Therm) Shift : " << i << ", Switch : " << j << ", Parity : Odd\n";
            }
            active_parity = odd;
            mpi::exchange::exchange_halos_cascade(field, geo, topo);
            plaquette_odd = mpi::heatbathcb::samples(field, geo, topo, hp, rng, active_parity);
            //plaquette.insert(plaquette.end(), std::make_move_iterator(plaquette_odd.begin()),
                             std::make_move_iterator(plaquette_odd.end()));
        }

        // Random shift
        mpi::shift::random_shift(field, geo, halo_shift, topo, rng);
    }

    //Sampling
    if (topo.rank == 0){
        std::cout << "Sampling : " << rp.N_shift*rp.N_switch_eo*2*rp.hp.N_samples << " <P> samples, " << rp.N_shift/rp.N_shift_topo << " Q samples\n";
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

        //Measure topo
        if (rp.topo and (i % rp.N_shift_topo == 0)) {
            tQE_current = mpi::observables::topo_charge_flowed(field, geo, flow, topo,
                                                               rp.N_steps_gf, rp.N_rk_steps);
            tQE_tot.insert(tQE_tot.end(), std::make_move_iterator(tQE_current.begin()),
                           std::make_move_iterator(tQE_current.end()));
        }
    }

    //===========================Output======================================

    // Flatten the vector
    if (topo.rank == 0) {
        // Write the output
        int precision_filename = 1;
        std::string filename = "HBQCB_" + std::to_string(L * n_core_dims) + "b" +
                               io::format_double(hp.beta, precision_filename) + "Ns" +
                               std::to_string(rp.N_shift) + "Nsw" + std::to_string(rp.N_switch_eo) +
                               "Np" + std::to_string(hp.N_samples) + "c" +
                               std::to_string(rp.cold_start) + "Nswp" +
                               std::to_string(hp.N_sweeps) + "Nh" + std::to_string(hp.N_hits) + "_plaquette";
        int precision = 10;
        io::save_double(plaquette, filename, precision);


        std::string filename_tQE = "HBQCB_" + std::to_string(L * n_core_dims) + "b" +
                               io::format_double(hp.beta, precision_filename) + "Ns" +
                               std::to_string(rp.N_shift) + "Nsw" + std::to_string(rp.N_switch_eo) +
                               "Np" + std::to_string(hp.N_samples) + "c" +
                               std::to_string(rp.cold_start) + "Nswp" +
                               std::to_string(hp.N_sweeps) + "Nh" + std::to_string(hp.N_hits) + "_topo";
        io::save_topo(tQE_tot, filename_tQE, precision);
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
