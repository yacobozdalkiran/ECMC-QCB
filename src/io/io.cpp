//
// Created by ozdalkiran-l on 1/14/26.
//

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
extern "C" {
#include <lime.h>
}
#include "io.h"

namespace fs = std::filesystem;

// Saves a vector of doubles in ../data/filename.txt
void io::save_double(const std::vector<double>& data, const std::string& filename, int precision) {
    // Create a data folder if doesn't exists
    fs::path dir("data");

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_plaquette.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::fixed << std::setprecision(precision);
    for (const double& x : data) {
        file << x << "\n";
    }
    file.close();
    std::cout << "Plaquette written in " << filepath << "\n";
}

void io::save_topo(const std::vector<double>& tQE, const std::string& filename, int precision) {
    // Create a data folder if doesn't exists
    fs::path dir("data");

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_topo.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::fixed << std::setprecision(precision);
    for (size_t i = 0; i < tQE.size(); i += 3) {
        file << tQE[i] << " " << tQE[i + 1] << " " << tQE[i + 2] << "\n";
    }
    file.close();
    std::cout << "Topology written in " << filepath << "\n";
};

void io::save_seed(std::mt19937_64& rng, const std::string& filename, mpi::MpiTopology& topo) {
    // Create a data folder if doesn't exists
    fs::path base_dir("data");
    fs::path run_dir =
        base_dir / (filename+"_seed");  // Utilise l'opérateur / pour gérer les slashs proprement

    try {
        // create_directories crée "data" PUIS "data/run_name" si nécessaire
        if (!fs::exists(run_dir)) {
            fs::create_directories(run_dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating directory structure " << run_dir << " : " << e.what()
                  << std::endl;
        return;
    }
    fs::path filepath = run_dir / (filename + "_seed" + std::to_string(topo.rank) + ".txt");

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << rng;
    file.close();
    if (topo.rank == 0) {
        std::cout << "Seed saved in " << filepath << "\n";
    }
};

void io::save_params(const RunParamsHbCB& rp, const std::string& filename) {
    // Create a data folder if doesn't exists
    fs::path dir("data");

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_params.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::boolalpha;

    file << "# Lattice params\n";
    file << "L_core=" << rp.L_core << "\n";
    file << "n_core_dims=" << rp.n_core_dims << "\n";
    file << "cold_start=" << rp.cold_start << "\n";
    file << "N_shift=" << rp.N_shift << "\n";
    file << "N_switch_eo=" << rp.N_switch_eo << "\n";
    file << "seed=" << rp.seed << "\n\n";

    file << "# Hb params\n";
    file << "beta = " << rp.hp.beta << "\n";
    file << "N_samples = " << rp.hp.N_samples << "\n";
    file << "N_sweeps = " << rp.hp.N_sweeps << "\n";
    file << "N_hits = " << rp.hp.N_hits << "\n\n";

    file << "# Run and topo params\n";
    file << "N_therm = " << rp.N_therm << "\n";
    file << "topo = " << rp.topo << "\n";
    file << "N_shift_topo = " << rp.N_shift_topo << "\n";
    file << "N_steps_gf = " << rp.N_steps_gf << "\n";
    file << "N_rk_steps = " << rp.N_rk_steps << "\n\n";
    file << "#########################################################\n\n";

    file.close();
    std::cout << "Parameters saved in " << filepath << "\n";
};

// Utilitary function to trim the spaces
std::string io::trim(const std::string& s) {
    size_t first = s.find_first_not_of(" \t");
    if (first == std::string::npos) return "";
    size_t last = s.find_last_not_of(" \t");
    return s.substr(first, (last - first + 1));
}

void io::load_params(const std::string& filename, RunParamsECB& rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Impossible d'ouvrir " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        // Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    // Lattice params
    if (config.count("L_core")) rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start")) rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed")) rp.seed = std::stoi(config["seed"]);
    if (config.count("N_switch_eo")) rp.N_switch_eo = std::stoi(config["N_switch_eo"]);
    if (config.count("N_shift")) rp.N_shift = std::stoi(config["N_shift"]);

    // ECMC params
    if (config.count("beta")) rp.ecmc_params.beta = std::stod(config["beta"]);
    if (config.count("N_samples")) rp.ecmc_params.N_samples = std::stoi(config["N_samples"]);
    if (config.count("param_theta_sample"))
        rp.ecmc_params.param_theta_sample = std::stod(config["param_theta_sample"]);
    if (config.count("param_theta_refresh"))
        rp.ecmc_params.param_theta_refresh = std::stod(config["param_theta_refresh"]);
    if (config.count("poisson")) rp.ecmc_params.poisson = (config["poisson"] == "true");
    if (config.count("epsilon_set")) rp.ecmc_params.epsilon_set = std::stod(config["epsilon_set"]);

    // Run and topo params
    if (config.count("N_therm")) rp.N_therm = std::stoi(config["N_therm"]);
    if (config.count("topo")) rp.topo = (config["topo"] == "true");
    if (config.count("N_shift_topo")) rp.N_shift_topo = std::stoi(config["N_shift_topo"]);
    if (config.count("N_steps_gf")) rp.N_steps_gf = std::stoi(config["N_steps_gf"]);
    if (config.count("N_rk_steps")) rp.N_rk_steps = std::stoi(config["N_rk_steps"]);
    //
    // Run name
    if (config.count("run_name")) rp.run_name = config["run_name"];
}

void io::load_params(const std::string& filename, RunParamsHbCB& rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open file " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        // Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    // Lattice params
    if (config.count("L_core")) rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start")) rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed")) rp.seed = std::stoi(config["seed"]);
    if (config.count("N_switch_eo")) rp.N_switch_eo = std::stoi(config["N_switch_eo"]);
    if (config.count("N_shift")) rp.N_shift = std::stoi(config["N_shift"]);

    // Hb params
    if (config.count("beta")) rp.hp.beta = std::stod(config["beta"]);
    if (config.count("N_samples")) rp.hp.N_samples = std::stoi(config["N_samples"]);
    if (config.count("N_sweeps")) rp.hp.N_sweeps = std::stoi(config["N_sweeps"]);
    if (config.count("N_hits")) rp.hp.N_hits = std::stoi(config["N_hits"]);

    // Run and topo params
    if (config.count("N_therm")) rp.N_therm = std::stoi(config["N_therm"]);
    if (config.count("topo")) rp.topo = (config["topo"] == "true");
    if (config.count("N_shift_topo")) rp.N_shift_topo = std::stoi(config["N_shift_topo"]);
    if (config.count("N_steps_gf")) rp.N_steps_gf = std::stoi(config["N_steps_gf"]);
    if (config.count("N_rk_steps")) rp.N_rk_steps = std::stoi(config["N_rk_steps"]);

    // Run name
    if (config.count("run_name")) rp.run_name = config["run_name"];
}

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

void print_time(long elapsed) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Elapsed time : " << elapsed << "s\n";
    std::cout << "==========================================" << std::endl;
}

std::string io::format_double(double val, int precision) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << val;
    return ss.str();
}
