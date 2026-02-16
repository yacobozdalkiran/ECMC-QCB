#ifndef GEOMETRYCB_H
#define GEOMETRYCB_H

#include <vector>

#include "types.h"

class GeometryCB {
public:
    int L;
    size_t V;
    size_t V_halo;

private:
    std::vector<size_t> neighbors;
    std::vector<bool> frozen;
    std::vector<std::pair<size_t, int>> links_staples;

public:
    explicit GeometryCB(int L_);

    // Index function for links
    [[nodiscard]] size_t index(int x, int y, int z, int t) const {
        return ((static_cast<size_t>(t) * L + z) * L + y) * L + x;
    }

    // Returns the index of the halo site corresponding to (c1,c2,c3) with order x,y,z,t
    [[nodiscard]] size_t index_halo_ecmc(int c1, int c2, int c3) const {
        return (static_cast<size_t>(c3) * L + c2) * L + c1;
    }

    // Index of a site taking into account halos
    [[nodiscard]] size_t index_w_halo(int x, int y, int z, int t) const {
        if (x == L) {
            return (V + index_halo_ecmc(y, z, t));
        }
        if (y == L) {
            return (V + V_halo + index_halo_ecmc(x, z, t));
        }
        if (z == L) {
            return (V + 2 * V_halo + index_halo_ecmc(x, y, t));
        }
        if (t == L) {
            return (V + 3 * V_halo + index_halo_ecmc(x, y, z));
        }
        if (x == -1) {
            return (V + 4 * V_halo + index_halo_ecmc(y, z, t));
        }
        if (y == -1) {
            return (V + 5 * V_halo + index_halo_ecmc(x, z, t));
        }
        if (z == -1) {
            return (V + 6 * V_halo + index_halo_ecmc(x, y, t));
        }
        if (t == -1) {
            return (V + 7 * V_halo + index_halo_ecmc(x, y, z));
        }
        return index(x, y, z, t);
    }

    // Index of a site in neighbor vector
    [[nodiscard]] static size_t index_neigh(size_t site, int mu, dir d) {
        return site * 8 + mu * 2 + d;
    }

    // Index of a link in is_frozen
    [[nodiscard]] static size_t index_frozen(size_t site, int mu) { return site * 4 + mu; }

    // Index of a link in links_staples
    [[nodiscard]] static size_t index_staples(size_t site, int mu, int i_staple, int i_link) {
        return site * 3 * 6 * 4 + mu * 3 * 6 + i_staple * 3 + i_link;
    }

    // Get a neighbor
    [[nodiscard]] size_t get_neigh(size_t site, int mu, dir d) const {
        return neighbors[index_neigh(site, mu, d)];
    }

    // Returns true if link is frozen
    [[nodiscard]] bool is_frozen(size_t site, int mu) const {
        return frozen[index_frozen(site, mu)];
    }

    // Returns the coords of the link of index i_link of the staple of index i_staple of <site,mu>
    [[nodiscard]] std::pair<size_t, int> get_link_staple(size_t site, int mu, int i_staple,
                                                         int i_link) const {
        return links_staples[index_staples(site, mu, i_staple, i_link)];
    }
};

#endif
