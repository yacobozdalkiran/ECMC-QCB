//
// Created by ozdalkiran-l on 1/13/26.
//

#include "ecmc_mpi_cb.h"

#include <iostream>

#include "../observables/observables_mpi.h"
#include "../su3/utils.h"

// Computes the list of the 6 staples around a gauge link
void mpi::ecmccb::compute_list_staples(const GaugeField& field, const GeometryCB& geo, size_t site,
                                       int mu, std::array<SU3, 6>& list_staple) {
    size_t index = 0;
    for (int nu = 0; nu < 4; nu++) {
        if (nu == mu) {
            continue;
        }
        size_t x = site;                              // x
        size_t xmu = geo.get_neigh(x, mu, up);        // x+mu
        size_t xnu = geo.get_neigh(x, nu, up);        // x+nu
        size_t xmunu = geo.get_neigh(xmu, nu, down);  // x+mu-nu
        size_t xmnu = geo.get_neigh(x, nu, down);     // x-nu
        auto U0 = field.view_link_const(xmu, nu);
        auto U1 = field.view_link_const(xnu, mu);
        auto U2 = field.view_link_const(x, nu);
        list_staple[index] = U0 * U1.adjoint() * U2.adjoint();
        auto V0 = field.view_link_const(xmunu, nu);
        auto V1 = field.view_link_const(xmnu, mu);
        auto V2 = field.view_link_const(xmnu, nu);
        list_staple[index + 1] = V0.adjoint() * V1.adjoint() * V2;
        index += 2;
    }
}

//Optimization of solve_reject
void mpi::ecmccb::solve_reject_fast(double A, double B, double& gamma, double& reject, int epsilon) {
    if (epsilon == -1) B = -B;

    double R = std::sqrt(A * A + B * B);
    double period = 2.0 * R;

    // Calcul du nombre de périodes rejetées (identique à l'original)
    double discarded_number = std::floor(gamma / period);
    gamma -= discarded_number * period;

    double phi = std::atan2(-A, B);
    if (phi < 0.0) phi += 2.0 * M_PI;

    double alpha;
    
    // Fusion des cas 1 (phi < pi/2) et 2 (phi > 3pi/2)
    if (phi < M_PI / 2.0 || phi > 3.0 * M_PI / 2.0) {
        double p1 = R - A;
        if (gamma > p1) {
            gamma -= p1;
            alpha = gamma / R - 1.0;
        } else {
            alpha = (gamma + A) / R;
        }
    } else {
        // Cas 3 (pi/2 <= phi <= 3pi/2)
        alpha = gamma / R - 1.0;
    }

    // Sécurité contre les erreurs d'arrondi flottant avant le std::asin
    alpha = std::max(-1.0, std::min(1.0, alpha));

    // Calcul de l'angle unique
    double theta = phi + std::asin(alpha);

    // Remplacement ultra-rapide de fmod(..., 2*M_PI)
    if (theta < 0.0) {
        theta += 2.0 * M_PI;
    } else if (theta >= 2.0 * M_PI) {
        theta -= 2.0 * M_PI;
    }

    reject = theta + 2.0 * M_PI * discarded_number;
}

// Solves the reject equation
void mpi::ecmccb::solve_reject(double A, double B, double& gamma, double& reject, int epsilon) {
    if (epsilon == -1) B = -B;
    double R = sqrt(A * A + B * B);
    double phi = atan2(-A, B);
    // Phi reduction
    if (phi < 0) phi += 2 * M_PI;
    double period = 0.0, p1 = 0.0, p2 = 0.0;
    // To store the 2 intervals on which the derivative of the action is positive
    std::array<double, 4> intervals = {0.0, 0.0, 2 * M_PI, 2 * M_PI};
    // Number of periods we discarded for the reject angle
    int discarded_number = 0;

    if (phi < M_PI / 2.0) {
        // cout << "cas 1"<< endl;
        intervals[1] = M_PI / 2.0 + phi;
        intervals[2] = 3 * M_PI / 2.0 + phi;
        p1 = R * (sin(intervals[1] - phi) - sin(intervals[0] - phi));
        p2 = R * (sin(intervals[3] - phi) - sin(intervals[2] - phi));
        if ((p1 < 0) && (p2 < 0)) std::cerr << "Périodes négatives !" << std::endl;
        period = p1 + p2;
        discarded_number = std::floor(gamma / period);
        gamma = gamma - std::floor(gamma / period) * period;
        if (gamma > p1) {
            gamma -= p1;
            double alpha = gamma / R + sin(intervals[2] - phi);
            double theta1 = fmod((phi + asin(alpha) + 2 * M_PI), 2 * M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2 * M_PI), 2 * M_PI);
            if ((theta1 < intervals[3]) && (theta1 > intervals[2])) {
                reject = theta1;
            } else {
                reject = theta2;
            }
        } else {
            double alpha = gamma / R + sin(intervals[0] - phi);
            double theta1 = fmod((phi + asin(alpha) + 2 * M_PI), 2 * M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2 * M_PI), 2 * M_PI);
            if ((theta1 < intervals[1]) && (theta1 > intervals[0])) {
                reject = theta1;
            } else {
                reject = theta2;
            }
        }
    }
    if (phi > 3 * M_PI / 2.0) {
        // cout << "cas 2" << endl;
        intervals[1] = -3 * M_PI / 2.0 + phi;
        intervals[2] = -M_PI / 2.0 + phi;
        // cout << "[" << intervals[0] << ", " << intervals[1] << "]" << endl;
        // cout << "[" << intervals[2] << ", " << intervals[3] << "]" << endl;
        p1 = R * (sin(intervals[1] - phi) - sin(intervals[0] - phi));
        p2 = R * (sin(intervals[3] - phi) - sin(intervals[2] - phi));
        if ((p1 < 0) && (p2 < 0)) std::cerr << "Périodes négatives !" << std::endl;
        period = p1 + p2;
        // cout << "contrib periodique = " << period << endl;
        discarded_number = std::floor(gamma / period);
        gamma = gamma - std::floor(gamma / period) * period;
        if (gamma > p1) {
            gamma -= p1;
            double alpha = gamma / R + sin(intervals[2] - phi);
            double theta1 = fmod((phi + asin(alpha) + 2 * M_PI), 2 * M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2 * M_PI), 2 * M_PI);
            if ((theta1 < intervals[3]) && (theta1 > intervals[2])) {
                reject = theta1;
            } else {
                reject = theta2;
            }
        } else {
            double alpha = gamma / R + sin(intervals[0] - phi);
            double theta1 = fmod((phi + asin(alpha) + 2 * M_PI), 2 * M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2 * M_PI), 2 * M_PI);
            if ((theta1 < intervals[1]) && (theta1 > intervals[0])) {
                reject = theta1;
            } else {
                reject = theta2;
            }
        }
    }
    if ((phi >= M_PI / 2.0) && (phi <= 3 * M_PI / 2.0)) {
        // cout << "cas 3" << endl;
        intervals[0] = -M_PI / 2.0 + phi;
        intervals[1] = M_PI / 2.0 + phi;
        period = R * (sin(intervals[1] - phi) - sin(intervals[0] - phi));
        if (period < 0) std::cerr << "Période négative !" << std::endl;
        // cout << "contrib periodique = " << period << endl;
        discarded_number = std::floor(gamma / period);
        gamma = gamma - std::floor(gamma / period) * period;
        double alpha = gamma / R + sin(intervals[0] - phi);
        double theta1 = fmod((phi + asin(alpha) + 2 * M_PI), 2 * M_PI);
        double theta2 = fmod((phi + M_PI - asin(alpha) + 2 * M_PI), 2 * M_PI);
        if ((theta1 < intervals[1]) && (theta1 > intervals[0])) {
            reject = theta1;
        } else {
            reject = theta2;
        }
    }
    reject += 2 * M_PI * discarded_number;
}

// Generates the 6 reject angles for a link
void mpi::ecmccb::compute_reject_angles(const GaugeField& field, size_t site, int mu,
                                        const std::array<SU3, 6>& list_staple, const SU3& R,
                                        int epsilon, const double& beta,
                                        std::array<double, 6>& reject_angles,
                                        std::mt19937_64& rng) {
    for (int i = 0; i < 6; i++) {
        std::uniform_real_distribution<double> unif01(0.0, 1.0);
        double gamma = -log(unif01(rng));
        SU3 P = R.adjoint() * field.view_link_const(site, mu) * list_staple[i] *
                R;  // * lambda_3 * Complex(0.0, 1.0);
        double A = P(0, 0).real() + P(1, 1).real();
        double B = -P(0, 0).imag() + P(1, 1).imag();
        A *= -(beta / 3.0);
        B *= -(beta / 3.0);
        solve_reject_fast(A, B, gamma, reject_angles[i], epsilon);
    }
}

// Selects an index between 0 and probas.size()-1 using the tower of probability method
size_t mpi::ecmccb::selectVariable(const std::array<double, 4>& probas, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double r = unif01(rng);
    double s = 0.0;
    for (size_t i = 0; i < 4; i++) {
        s += probas[i];
        if (s > r) {
            return i;
        }
    }
    std::cerr << "SelectVariable Error" << std::endl;
    return -1;
}

// Returns a new link, direction and R matrix for a lift
std::pair<std::pair<size_t, int>, int> mpi::ecmccb::lift_improved(
    const GaugeField& field, const GeometryCB& geo, size_t site, int mu, int j, SU3& R,
    const SU3& lambda_3, const std::vector<SU3>& set, std::mt19937_64& rng) {
    std::array<std::pair<size_t, int>, 4>
        links_plaquette_j;  // We add the current link to get the plaquette
    links_plaquette_j[0] = std::make_pair(site, mu);
    links_plaquette_j[1] = geo.get_link_staple(site, mu, j, 0);
    links_plaquette_j[2] = geo.get_link_staple(site, mu, j, 1);
    links_plaquette_j[3] = geo.get_link_staple(site, mu, j, 2);

    SU3 U0 = field.view_link_const(site, mu);
    SU3 U1 = field.view_link_const(links_plaquette_j[1].first, links_plaquette_j[1].second);
    SU3 U2 = field.view_link_const(links_plaquette_j[2].first, links_plaquette_j[2].second);
    SU3 U3 = field.view_link_const(links_plaquette_j[3].first, links_plaquette_j[3].second);
    std::array<double, 4> probas{};
    std::array<double, 4> abs_dS{};
    double sum = 0.0;
    std::array<int, 4> sign_dS{};
    std::array<SU3, 4> P{};

    if (j % 2 == 0) {  // Forward plaquette
        P[0] = U0 * U1 * U2.adjoint() * U3.adjoint();
        if (!geo.is_frozen(links_plaquette_j[1].first, links_plaquette_j[1].second))
            P[1] = U1 * U2.adjoint() * U3.adjoint() * U0;
        if (!geo.is_frozen(links_plaquette_j[2].first, links_plaquette_j[2].second))
            P[2] = U2 * U1.adjoint() * U0.adjoint() * U3;
        if (!geo.is_frozen(links_plaquette_j[3].first, links_plaquette_j[3].second))
            P[3] = U3 * U2 * U1.adjoint() * U0.adjoint();
    } else {  // Backward plaquette
        P[0] = U0 * U1.adjoint() * U2.adjoint() * U3;
        if (!geo.is_frozen(links_plaquette_j[1].first, links_plaquette_j[1].second))
            P[1] = U1 * U0.adjoint() * U3.adjoint() * U2;
        if (!geo.is_frozen(links_plaquette_j[2].first, links_plaquette_j[2].second))
            P[2] = U2 * U1 * U0.adjoint() * U3.adjoint();
        if (!geo.is_frozen(links_plaquette_j[3].first, links_plaquette_j[3].second))
            P[3] = U3 * U0 * U1.adjoint() * U2.adjoint();
    }
    for (size_t i = 0; i < 4; i++) {
        if (!geo.is_frozen(links_plaquette_j[i].first, links_plaquette_j[i].second))
            probas[i] = -(Complex(0.0, 1.0) * lambda_3 * R.adjoint() * P[i] * R).trace().real();
        else
            probas[i] = 0.0;
        sign_dS[i] = dsign(probas[i]);
        probas[i] = abs(probas[i]);
        abs_dS[i] = probas[i];
        sum += probas[i];
    }
    for (size_t i = 0; i < 4; i++) {
        probas[i] /= sum;
    }
    size_t index_lift = selectVariable(probas, rng);

    // Change R
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    std::uniform_int_distribution<size_t> set_index(0, set.size() - 1);
    SU3 R_new;
    size_t i_set = set_index(rng);
    R_new = set[i_set] * R;
    double dS_j_R =
        (Complex(0.0, -1.0) * lambda_3 * R.adjoint() * P[index_lift] * R).trace().real();
    double dS_j_R_new =
        (Complex(0.0, -1.0) * lambda_3 * R_new.adjoint() * P[index_lift] * R_new).trace().real();
    double new_epsilon = -sign_dS[index_lift];
    if (abs(dS_j_R) < abs(dS_j_R_new)) {
        R = R_new;
        new_epsilon = -dsign(dS_j_R_new);
    } else {
        double r = unif01(rng);
        if (r < abs(dS_j_R_new) / abs(dS_j_R)) {
            R = R_new;
            new_epsilon = -dsign(dS_j_R_new);
        }
    }
    return make_pair(links_plaquette_j[index_lift], new_epsilon);
}

// Updates the gauge field with XY embedding
void mpi::ecmccb::update(GaugeField& field, size_t site, int mu, double theta, int epsilon,
                         const SU3& R) {
    SU3 Uold = field.view_link_const(site, mu);
    field.view_link(site, mu) = R * el_3(epsilon * theta) * R.adjoint() * Uold;
    // projection_su3(links, site, mu);
}

// Returns a random non frozen site
size_t mpi::ecmccb::random_site(const GeometryCB& geo, std::mt19937_64& rng) {
    int L = geo.L_int;
    std::uniform_int_distribution random_coord(1, L - 1);
    int x = random_coord(rng);
    int y = random_coord(rng);
    int z = random_coord(rng);
    int t = random_coord(rng);
    return geo.index(x, y, z, t);
}

// Generates samples of global mean plaquette using ECMC with frozen BC
std::vector<double> mpi::ecmccb::samples_improved(GaugeField& field, const GeometryCB& geo,
                                                  const ECMCParams& params, std::mt19937_64& rng,
                                                  mpi::MpiTopology& topo, parity active_parity) {
    double beta = params.beta;
    int N_samples = params.N_samples;
    double param_theta_sample = params.param_theta_sample;
    double param_theta_refresh = params.param_theta_refresh;
    bool poisson = params.poisson;
    double epsilon_set = params.epsilon_set;

    // if (param_theta_sample<param_theta_refresh) {
    //     std::cerr << "Wrong args value, must have param_theta_sample>param_theta_refresh \n";
    // }
    // Set de matrices pour refresh R
    int N_set = 100;
    size_t lift_counter = 0;
    std::vector<SU3> set(N_set + 1);
    ecmc_set(epsilon_set, set, rng);

    // Variables aléatoires
    std::uniform_int_distribution<int> random_dir(0, 3);
    std::uniform_int_distribution<int> random_eps(0, 1);
    std::exponential_distribution<double> random_theta_sample(1.0 / param_theta_sample);
    std::exponential_distribution<double> random_theta_refresh(1.0 / param_theta_refresh);

    // Matrice lambda_3 de Gell-Mann
    SU3 lambda_3;
    lambda_3 << Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0),
        Complex(-1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0),
        Complex(0.0, 0.0);

    // Initialisation aléatoire de la position de la chaîne
    size_t site_current = random_site(geo, rng);
    int mu_current = random_dir(rng);
    int epsilon_current = 2 * random_eps(rng) - 1;

    // Initialisation aléatoire des theta limites pour sample et refresh
    double theta_sample{};
    double theta_refresh{};
    if (poisson) {
        theta_sample = random_theta_sample(rng);
        theta_refresh = random_theta_refresh(rng);
    } else {
        theta_sample = param_theta_sample;
        theta_refresh = param_theta_refresh;
    }

    // Initialisation des angles totaux parcourus à 0.0
    double theta_parcouru_sample = 0.0;
    double theta_parcouru_refresh = 0.0;

    // Angle d'update
    double theta_update = 0.0;

    // Arrays utilisés à chaque étape de la chaîne (évite de les initialiser des milliers de fois)
    std::array<double, 6> reject_angles = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<SU3, 6> list_staple;

    SU3 R = random_su3(rng);

    int samples = 0;
    std::array<double, 2> deltas = {0.0, 0.0};
    // size_t event_counter = 0;
    std::vector<double> meas_plaquette{};
    if (topo.p == active_parity) {
        while (samples < N_samples) {
            if (lift_counter % N_set == 0) {
                ecmc_set(epsilon_set, set, rng);
            }
            compute_list_staples(field, geo, site_current, mu_current, list_staple);
            compute_reject_angles(field, site_current, mu_current, list_staple, R, epsilon_current,
                                  beta, reject_angles, rng);
            auto it = std::min_element(reject_angles.begin(), reject_angles.end());
            auto j = static_cast<int>(
                std::distance(reject_angles.begin(), it));  // theta_reject = reject_angles[j]
            // cout << "Angle reject : " << reject_angles[j] << endl;
            deltas[0] = theta_sample - reject_angles[j] - theta_parcouru_sample;
            deltas[1] = theta_refresh - reject_angles[j] - theta_parcouru_refresh;
            if (deltas[1] > 0) {
                // lifts++;
            }
            // proposed++;

            auto it_deltas = std::min_element(deltas.begin(), deltas.end());
            auto F = static_cast<int>(std::distance(deltas.begin(), it_deltas));

            if ((deltas[0] < 0) && (deltas[1] < 0)) {
                if (F == 0) {
                    // On sample
                    theta_update = theta_sample - theta_parcouru_sample;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    double plaq = mpi::observables::mean_plaquette_global(field, geo, topo);
                    if (topo.rank == 0) {
                        std::cout << "Sample " << samples << ", ";
                        std::cout << "<P> = " << plaq << "\n";
                    }
                    meas_plaquette.emplace_back(plaq);
                    samples++;
                    theta_parcouru_sample = 0;
                    if (poisson)
                        theta_sample =
                            random_theta_sample(rng);  // On retire un nouveau theta_sample
                    theta_parcouru_refresh += theta_update;
                    // On update jusqu'au refresh
                    theta_update = theta_refresh - theta_parcouru_refresh;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    theta_parcouru_sample += theta_update;
                    theta_parcouru_refresh = 0;
                    if (poisson)
                        theta_refresh = random_theta_refresh(rng);  // On retire un nouveau theta
                                                                    // refresh
                    // On refresh
                    // event_counter++;
                    site_current = random_site(geo, rng);
                    mu_current = random_dir(rng);
                    epsilon_current = 2 * random_eps(rng) - 1;
                    R = random_su3(rng);
                }
                if (F == 1) {
                    // On update jusqu'au refresh
                    theta_update = theta_refresh - theta_parcouru_refresh;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    theta_parcouru_sample += theta_update;
                    theta_parcouru_refresh = 0;
                    if (poisson)
                        theta_refresh = random_theta_refresh(rng);  // On retire un nouveau
                                                                    // theta_refresh
                    // On refresh
                    // event_counter++;
                    site_current = random_site(geo, rng);
                    mu_current = random_dir(rng);
                    epsilon_current = 2 * random_eps(rng) - 1;
                    R = random_su3(rng);
                }
            } else if (deltas[F] < 0) {
                if (F == 0) {
                    // On update jusqu'a theta_sample
                    theta_update = theta_sample - theta_parcouru_sample;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    // On sample
                    double plaq = mpi::observables::mean_plaquette_global(field, geo, topo);
                    if (topo.rank == 0) {
                        std::cout << "Sample " << samples << ", ";
                        std::cout << "<P> = " << plaq << "\n";
                    }
                    meas_plaquette.emplace_back(plaq);
                    samples++;
                    theta_parcouru_sample = 0;
                    if (poisson)
                        theta_sample =
                            random_theta_sample(rng);  // On retire un nouveau theta_sample
                    theta_parcouru_refresh += theta_update;
                    // On finit l'update et on lift
                    theta_update = -deltas[F];
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    theta_parcouru_sample += theta_update;
                    theta_parcouru_refresh += theta_update;
                    // On lifte
                    // event_counter++;
                    auto l = lift_improved(field, geo, site_current, mu_current, j, R, lambda_3,
                                           set, rng);
                    lift_counter++;
                    site_current = l.first.first;
                    mu_current = l.first.second;
                    epsilon_current = l.second;
                }
                if (F == 1) {
                    // On update jusqu'à theta_refresh
                    theta_update = theta_refresh - theta_parcouru_refresh;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    theta_parcouru_sample += theta_update;
                    theta_parcouru_refresh = 0;
                    if (poisson)
                        theta_refresh = random_theta_refresh(rng);  // On retire un nouveau
                                                                    // theta_refresh
                    // On refresh
                    // event_counter++;
                    site_current = random_site(geo, rng);
                    mu_current = random_dir(rng);
                    epsilon_current = 2 * random_eps(rng) - 1;
                    R = random_su3(rng);
                }
            } else {
                // On update
                theta_update = reject_angles[j];
                update(field, site_current, mu_current, theta_update, epsilon_current, R);
                theta_parcouru_sample += theta_update;
                theta_parcouru_refresh += theta_update;
                // On lift
                // event_counter++;
                auto l =
                    lift_improved(field, geo, site_current, mu_current, j, R, lambda_3, set, rng);
                lift_counter++;
                site_current = l.first.first;
                mu_current = l.first.second;
                epsilon_current = l.second;
            }
        }
    } else {
        while (samples < N_samples) {
            double plaq = mpi::observables::mean_plaquette_global(field, geo, topo);
            if (topo.rank == 0) {
                std::cout << "Sample " << samples << ", ";
                std::cout << "<P> = " << plaq << "\n";
            }
            meas_plaquette.emplace_back(plaq);
            samples++;
        }
    }

    return meas_plaquette;
}

// Generates 1 of using ECMC with frozen BC
void mpi::ecmccb::sample(GaugeField& field, const GeometryCB& geo, const ECMCParams& params,
                         std::mt19937_64& rng, mpi::MpiTopology& topo, parity active_parity) {
    double beta = params.beta;
    double param_theta_sample = params.param_theta_sample;
    double param_theta_refresh = params.param_theta_refresh;
    bool poisson = params.poisson;
    double epsilon_set = params.epsilon_set;

    // if (param_theta_sample<param_theta_refresh) {
    //     std::cerr << "Wrong args value, must have param_theta_sample>param_theta_refresh \n";
    // }
    // Set de matrices pour refresh R
    int N_set = 100;
    size_t lift_counter = 0;
    std::vector<SU3> set(N_set + 1);
    ecmc_set(epsilon_set, set, rng);

    // Variables aléatoires
    std::uniform_int_distribution<int> random_dir(0, 3);
    std::uniform_int_distribution<int> random_eps(0, 1);
    std::exponential_distribution<double> random_theta_sample(1.0 / param_theta_sample);
    std::exponential_distribution<double> random_theta_refresh(1.0 / param_theta_refresh);

    // Matrice lambda_3 de Gell-Mann
    SU3 lambda_3;
    lambda_3 << Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0),
        Complex(-1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0),
        Complex(0.0, 0.0);

    // Initialisation aléatoire de la position de la chaîne
    size_t site_current = random_site(geo, rng);
    int mu_current = random_dir(rng);
    int epsilon_current = 2 * random_eps(rng) - 1;

    // Initialisation aléatoire des theta limites pour sample et refresh
    double theta_sample{};
    double theta_refresh{};
    if (poisson) {
        theta_sample = random_theta_sample(rng);
        theta_refresh = random_theta_refresh(rng);
    } else {
        theta_sample = param_theta_sample;
        theta_refresh = param_theta_refresh;
    }

    // Initialisation des angles totaux parcourus à 0.0
    double theta_parcouru_sample = 0.0;
    double theta_parcouru_refresh = 0.0;

    // Angle d'update
    double theta_update = 0.0;

    // Arrays utilisés à chaque étape de la chaîne (évite de les initialiser des milliers de fois)
    std::array<double, 6> reject_angles = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<SU3, 6> list_staple;

    SU3 R = random_su3(rng);

    std::array<double, 2> deltas = {0.0, 0.0};
    // size_t event_counter = 0;
    std::vector<double> meas_plaquette{};
    if (topo.p == active_parity) {
        while (1) {
            if (lift_counter % N_set == 0) {
                ecmc_set(epsilon_set, set, rng);
            }
            compute_list_staples(field, geo, site_current, mu_current, list_staple);
            compute_reject_angles(field, site_current, mu_current, list_staple, R, epsilon_current,
                                  beta, reject_angles, rng);
            auto it = std::min_element(reject_angles.begin(), reject_angles.end());
            auto j = static_cast<int>(
                std::distance(reject_angles.begin(), it));  // theta_reject = reject_angles[j]
            // cout << "Angle reject : " << reject_angles[j] << endl;
            deltas[0] = theta_sample - reject_angles[j] - theta_parcouru_sample;
            deltas[1] = theta_refresh - reject_angles[j] - theta_parcouru_refresh;
            if (deltas[1] > 0) {
                // lifts++;
            }
            // proposed++;

            auto it_deltas = std::min_element(deltas.begin(), deltas.end());
            auto F = static_cast<int>(std::distance(deltas.begin(), it_deltas));

            if ((deltas[0] < 0) && (deltas[1] < 0)) {
                if (F == 0) {
                    // On update jusqu'au sample
                    theta_update = theta_sample - theta_parcouru_sample;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    return;
                }
                if (F == 1) {
                    // On update jusqu'au refresh
                    theta_update = theta_refresh - theta_parcouru_refresh;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    theta_parcouru_sample += theta_update;
                    theta_parcouru_refresh = 0;
                    if (poisson)
                        theta_refresh = random_theta_refresh(rng);  // On retire un nouveau
                                                                    // theta_refresh
                    // On refresh
                    // event_counter++;
                    site_current = random_site(geo, rng);
                    mu_current = random_dir(rng);
                    epsilon_current = 2 * random_eps(rng) - 1;
                    R = random_su3(rng);
                }
            } else if (deltas[F] < 0) {
                if (F == 0) {
                    // On update jusqu'a theta_sample
                    theta_update = theta_sample - theta_parcouru_sample;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    return;
                }
                if (F == 1) {
                    // On update jusqu'à theta_refresh
                    theta_update = theta_refresh - theta_parcouru_refresh;
                    update(field, site_current, mu_current, theta_update, epsilon_current, R);
                    theta_parcouru_sample += theta_update;
                    theta_parcouru_refresh = 0;
                    if (poisson)
                        theta_refresh = random_theta_refresh(rng);  // On retire un nouveau
                                                                    // theta_refresh
                    // On refresh
                    // event_counter++;
                    site_current = random_site(geo, rng);
                    mu_current = random_dir(rng);
                    epsilon_current = 2 * random_eps(rng) - 1;
                    R = random_su3(rng);
                }
            } else {
                // On update
                theta_update = reject_angles[j];
                update(field, site_current, mu_current, theta_update, epsilon_current, R);
                theta_parcouru_sample += theta_update;
                theta_parcouru_refresh += theta_update;
                // On lift
                // event_counter++;
                auto l =
                    lift_improved(field, geo, site_current, mu_current, j, R, lambda_3, set, rng);
                lift_counter++;
                site_current = l.first.first;
                mu_current = l.first.second;
                epsilon_current = l.second;
            }
        }
    }
}
