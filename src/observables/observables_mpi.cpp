//
// Created by ozdalkiran-l on 1/8/26.
//

#include "observables_mpi.h"
#include "../mpi/HalosExchange.h"

//Computation of mean plaquette with halos embedded in field (needs field halos exchange first)
double mpi::observables::mean_plaquette_local(const GaugeField& field, const GeometryCB& geo) {
    double sum = 0.0;
    SU3 U1, U2, U3, U4;
    for (int t = 1; t <= geo.L_int; t++) {
        for (int z = 1; z <= geo.L_int; z++) {
            for (int y = 1; y <= geo.L_int; y++) {
                for (int x = 1; x <= geo.L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        for (int nu = mu + 1; nu < 4; nu++) {
                            U1 = field.view_link_const(site, mu);
                            U2 = field.view_link_const(geo.get_neigh(site, mu, up), nu);
                            U3 = field.view_link_const(geo.get_neigh(site, nu, up), mu).adjoint();
                            U4 = field.view_link_const(site, nu).adjoint();
                            sum += (U1 * U2 * U3 * U4).trace().real() / 3.0;
                        }
                    }
                }
            }
        }
    }
    return sum;
}

//Computation of global mean plaquette with halos embedded in field (needs field halos exchanges first)
double mpi::observables::mean_plaquette_global(GaugeField &field, const GeometryCB &geo, MpiTopology &topo){
    double local_mean_plaquette = mean_plaquette_local(field, geo);
    double global_mean_plaquette = 0.0;
    MPI_Allreduce(&local_mean_plaquette, &global_mean_plaquette, 1, MPI_DOUBLE, MPI_SUM,
                  topo.cart_comm);
    global_mean_plaquette /= 6.0 * geo.V_int * topo.size;
    return global_mean_plaquette;
}
