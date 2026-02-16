//
// Created by ozdalkiran-l on 1/8/26.
//

#ifndef INC_4D_MPI_OBSERVABLES_H
#define INC_4D_MPI_OBSERVABLES_H

#include "../gauge/GaugeField.h"
#include "../mpi/MpiTopology.h"


namespace mpi::observables{
    double mean_plaquette_local(const GaugeField &field, const GeometryCB &geo);
    double mean_plaquette_global(GaugeField &field, const GeometryCB &geo, MpiTopology &topo);
}

#endif //INC_4D_MPI_OBSERVABLES_H
