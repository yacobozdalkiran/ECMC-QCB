//
//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_IO_H
#define ECMC_MPI_IO_H

#include <random>
#include <vector>

#include "../mpi/MpiTopology.h"
#include "params.h"

namespace io {
// Output
void save_double(const std::vector<double>& data, const std::string& filename, int precision);
void save_topo(const std::vector<double>& tQE, const std::string& filename, int precision);
void save_seed(std::mt19937_64& rng, const std::string& filename, mpi::MpiTopology& topo);
void save_params(const RunParamsHbCB& rp, const std::string& filename);
void save_params(const RunParamsECB& rp, const std::string& filename);
// Input
std::string trim(const std::string& s);
void load_params(const std::string& filename, RunParamsECB& rp);
void load_params(const std::string& filename, RunParamsHbCB& rp);

// Utilitaries
std::string format_double(double val, int precision);

//Load parameters from input file
bool read_params(RunParamsHbCB& params, int rank, const std::string& input);
bool read_params(RunParamsECB& params, int rank, const std::string& input);
}  // namespace io

// Printing
void print_parameters(const RunParamsHbCB& rp, const mpi::MpiTopology& topo);
void print_parameters(const RunParamsECB& rp, const mpi::MpiTopology& topo);
void print_time(long elapsed);

#endif  // ECMC_MPI_IO_H
