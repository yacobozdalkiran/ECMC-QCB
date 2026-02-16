#include "GeometryCB.h"

#include <cstdint>

GeometryCB::GeometryCB(int L_) {
    L = L_;
    V = L * L * L * L;
    V_halo = L * L * L;
    neighbors.resize((V + 8 * V_halo) * 8, SIZE_MAX);
    // On étend la boucle de -1 à L pour couvrir la "coquille" des halos
    for (int t = -1; t <= L; t++) {
        for (int z = -1; z <= L; z++) {
            for (int y = -1; y <= L; y++) {
                for (int x = -1; x <= L; x++) {
                    // On ne calcule les voisins que pour les sites qui existent dans notre
                    // indexation (index_w_halo ne gère pas les coins comme x=-1 ET y=-1, seulement
                    // une face à la fois)
                    int out_count = (x == -1 || x == L) + (y == -1 || y == L) +
                                    (z == -1 || z == L) + (t == -1 || t == L);

                    if (out_count > 1) continue;  // On ignore les coins/arêtes (pas nécessaires)

                    size_t site_idx = index_w_halo(x, y, z, t);

                    // Direction X (mu=0)
                    neighbors[index_neigh(site_idx, 0, up)] = index_w_halo(x + 1, y, z, t);
                    neighbors[index_neigh(site_idx, 0, down)] = index_w_halo(x - 1, y, z, t);
                    // Direction Y (mu=1)
                    neighbors[index_neigh(site_idx, 1, up)] = index_w_halo(x, y + 1, z, t);
                    neighbors[index_neigh(site_idx, 1, down)] = index_w_halo(x, y - 1, z, t);
                    // Direction Z (mu=2)
                    neighbors[index_neigh(site_idx, 2, up)] = index_w_halo(x, y, z + 1, t);
                    neighbors[index_neigh(site_idx, 2, down)] = index_w_halo(x, y, z - 1, t);
                    // Direction T (mu=3)
                    neighbors[index_neigh(site_idx, 3, up)] = index_w_halo(x, y, z, t + 1);
                    neighbors[index_neigh(site_idx, 3, down)] = index_w_halo(x, y, z, t - 1);
                }
            }
        }
    }

    frozen.resize((V + 8 * V_halo) * 4, false);
    // Frozen links are those that step out of the lattice core or belong to halos
    for (int t = -1; t <= L; t++) {
        for (int z = -1; z <= L; z++) {
            for (int y = -1; y <= L; y++) {
                for (int x = -1; x <= L; x++) {
                    bool link_is_frozen = false;
                    if (x == -1 or y == -1 or z == -1 or t == -1) link_is_frozen = true;
                    if (x == L or y == L or z == L or t == L) link_is_frozen = true;
                    for (int mu = 0; mu < 4; mu++) {
                        if (x == L - 1 and mu == 0) link_is_frozen = true;
                        if (y == L - 1 and mu == 1) link_is_frozen = true;
                        if (z == L - 1 and mu == 2) link_is_frozen = true;
                        if (t == L - 1 and mu == 3) link_is_frozen = true;
                        size_t i = index_w_halo(x, y, z, t);
                        frozen[index_frozen(i, mu)] = link_is_frozen;
                    }
                }
            }
        }
    }

    links_staples.resize(V * 4 * 6 * 3, std::make_pair(SIZE_MAX, -1));
    // 4 links per site, 6 staples per link, 3 links per staple
    for (int t = 0; t < L; t++) {
        for (int z = 0; z < L; z++) {
            for (int y = 0; y < L; y++) {
                for (int x = 0; x < L; x++) {
                    size_t site = index(x, y, z, t);  // x
                    for (int mu = 0; mu < 4; mu++) {
                        if (!is_frozen(site, mu)) {
                            int j = 0;
                            for (int nu = 0; nu < 4; nu++) {
                                if (nu == mu) continue;

                                size_t xmu = get_neigh(site, mu, up);     // x+mu
                                size_t xnu = get_neigh(site, nu, up);     // x+nu
                                size_t xmunu = get_neigh(xmu, nu, down);  // x+mu-nu
                                size_t xmnu = get_neigh(site, nu, down);  // x-nu

                                links_staples[index_staples(site, mu, j, 0)] = {xmu, nu};
                                links_staples[index_staples(site, mu, j, 1)] = {xnu, mu};
                                links_staples[index_staples(site, mu, j, 2)] = {site, nu};
                                links_staples[index_staples(site, mu, j + 1, 0)] = {xmunu, nu};
                                links_staples[index_staples(site, mu, j + 1, 1)] = {xmnu, mu};
                                links_staples[index_staples(site, mu, j + 1, 2)] = {xmnu, nu};

                                j += 2;
                            }
                        }
                    }
                }
            }
        }
    }
}
