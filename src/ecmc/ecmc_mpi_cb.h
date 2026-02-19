//
// Created by ozdalkiran-l on 1/13/26.
//

#ifndef INC_4D_MPI_ECMC_MPI_H
#define INC_4D_MPI_ECMC_MPI_H

#include <random>

#include "../gauge/GaugeField.h"
#include "../io/params.h"
#include "../mpi/MpiTopology.h"

namespace mpi::ecmccb {
void compute_list_staples(const GaugeField& field, const GeometryCB& geo, size_t site, int mu,
                          std::array<SU3, 6>& list_staple);
void solve_reject_fast(double A, double B, double& gamma, double& reject, int epsilon);
void compute_reject_angles(const GaugeField& field, size_t site, int mu,
                           const std::array<SU3, 6>& list_staple, const SU3& R, int epsilon,
                           const double& beta, std::array<double, 6>& reject_angles,
                           std::mt19937_64& rng);
size_t selectVariable(const std::array<double, 4>& probas, std::mt19937_64& rng);
double compute_ds(const SU3& Pi, const SU3& R_mat);
std::pair<std::pair<size_t, int>, int> lift_improved_fast(const GaugeField& field,
                                                          const GeometryCB& geo, size_t site,
                                                          int mu, int j, SU3& R,
                                                          const std::vector<SU3>& set,
                                                          std::mt19937_64& rng);
void update(GaugeField& field, size_t site, int mu, double theta, int epsilon, const SU3& R);
size_t random_site(const GeometryCB& geo, std::mt19937_64& rng);
void sample(GaugeField& field, const GeometryCB& geo, const ECMCParams& params,
            std::mt19937_64& rng, mpi::MpiTopology& topo, parity active_parity);
}  // namespace mpi::ecmccb

#endif  // INC_4D_MPI_ECMC_MPI_H
