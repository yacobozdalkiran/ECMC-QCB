//
//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_IO_H
#define ECMC_MPI_IO_H

#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

#include "params.h"

namespace io {
    //Output
    void save_double(const std::vector<double> &data, const std::string &filename, int precision);
    void save_double_params(const std::vector<double> &data, const RunParams &params, const std::string &filename, int precision);

    //Input
    std::string trim(const std::string& s);
    void load_params(const std::string& filename, RunParams& rp);
    void load_params(const std::string& filename, RunParamsCB& rp);
    void load_params(const std::string &filename, RunParamsSC& rp);
    void load_params(const std::string &filename, RunParamsMetro& rp);
    void load_params(const std::string &filename, RunParamsHb& rp);
    void load_params(const std::string &filename, RunParamsHbMPI& rp);
    void load_params(const std::string &filename, RunParamsECB& rp);
    void load_params(const std::string &filename, RunParamsHbCB& rp);

    //Utilitaries
    inline std::string format_double(double val, int precision) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(precision) << val;
        return ss.str();
    }

    //ILDG
    namespace ildg {
        std::string generate_ildg_xml(int lx, int ly, int lz, int lt, int precision = 64);
    }
}

//Printing

//Prints the elapsed time
inline void print_time(long elapsed) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Elapsed time : " << elapsed << "s\n";
    std::cout << "==========================================" << std::endl;
}


#endif //ECMC_MPI_IO_H
