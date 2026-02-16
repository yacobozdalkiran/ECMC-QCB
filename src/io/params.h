//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_PARAMS_H
#define ECMC_MPI_PARAMS_H


struct ECMCParams {
    double beta = 6.0;
    int N_samples = 10;
    double param_theta_sample = 100;
    double param_theta_refresh = 30;
    bool poisson = false;
    double epsilon_set = 0.15;
};

struct MetroParams {
    double beta = 6.0;
    double epsilon = 0.15;
    int N_set = 50; //Number of sweeps between set change
    int N_samples = 10; //Number of measures (samples)
    int N_sweeps_meas = 1; //Number of sweeps between measures
    int N_hits = 1;
    int N_burnin = 0;
};

struct HbParams {
    double beta = 6.0;
    int N_samples = 10;
    int N_hits = 1;
    int N_sweeps = 1;
};

struct RunParams {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int L_shift = 2;
    int n_shift = 1; //number of shifts in all the directions
    bool stype_pos=true; //Are the shifts positives ?
    int seed = 123;
    ECMCParams ecmc_params{};
};


struct RunParamsECB {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int seed = 123;
    int N_switch_eo = 3; //Number of switch even/odd per shift
    int N_shift = 2; //Number of shifts
    ECMCParams ecmc_params{}; //Params of the ECMC for each even/odd update
};

struct RunParamsCB {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int n_shift = 1; //number of shifts in all the directions
    bool stype_pos=true; //Are the shifts positives ?
    int seed = 123;
    ECMCParams ecmc_params{}; //Params for each color update (total samples nb = N_samples*2*N_shift)
};

struct RunParamsSC {
    int L=4;
    int T=4;
    bool cold_start=true;
    int seed = 123;
    ECMCParams ecmc_params{};
};

struct RunParamsMetro {
    int L = 4;
    int T = 4;
    bool cold_start=true;
    int seed = 123;
    MetroParams mp{};
};

struct RunParamsHb {
    int L=4;
    int T=4;
    bool cold_start=true;
    int seed = 123;
    HbParams hp{};
};

struct RunParamsHbMPI {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int L_shift = 2;
    int n_shift = 1; //number of shifts in all the directions
    bool stype_pos=true; //Are the shifts positives ?
    int seed = 123;
    HbParams hp{};
};


struct RunParamsHbCB {
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;
    int N_switch_eo = 3; //Number of switch even/odd per shift
    int N_shift = 2; //Number of shifts
    int seed = 123;
    HbParams hp{}; // Hb params for each even/odd update
};

#endif //ECMC_MPI_PARAMS_H
