#include <iostream>
#include <random>
#include <chrono>
#include <vector>

#include "../ecmc/ecmc_mpi_cb.h"
#include "../io/io.h"


int main() {
    std::random_device rd;
    std::mt19937_64 rng(rd());

    std::uniform_real_distribution<double> dist_coeffs(-2.0, 2.0);
    std::uniform_real_distribution<double> dist_gamma(0.01, 10.0);
    std::uniform_int_distribution<int> dist_eps(0, 1);

    int N = 10000000; // Augmenté pour plus de stabilité statistique
    
    // On pré-génère les données pour ne pas inclure le coût de la RNG dans le bench
    struct Params { double A, B, gamma; int eps; };
    std::vector<Params> data(N);
    for (int i = 0; i < N; ++i) {
        data[i] = {dist_coeffs(rng), dist_coeffs(rng), dist_gamma(rng), 2 * dist_eps(rng) - 1};
    }

    double dummy_reject = 0; // Pour éviter que le compilateur n'optimise et supprime l'appel

    // --- Benchmark Slow ---
    auto start_slow = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++) {
        mpi::ecmccb::solve_reject(data[i].A, data[i].B, data[i].gamma, dummy_reject, data[i].eps);
    }
    auto end_slow = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_slow = end_slow - start_slow;

    // --- Benchmark Fast ---
    auto start_fast = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++) {
        mpi::ecmccb::solve_reject_fast(data[i].A, data[i].B, data[i].gamma, dummy_reject, data[i].eps);
    }
    auto end_fast = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_fast = end_fast - start_fast;

    // --- Résultats ---
    std::cout << "Résultats sur " << N << " itérations :" << std::endl;
    std::cout << "Temps total Slow : " << duration_slow.count() << " ms" << std::endl;
    std::cout << "Temps total Fast : " << duration_fast.count() << " ms" << std::endl;
    
    double speedup = duration_slow.count() / duration_fast.count();
    std::cout << "Speedup : x" << speedup << std::endl;

    return 0;
}
