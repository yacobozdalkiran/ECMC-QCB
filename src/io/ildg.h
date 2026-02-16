#include <arpa/inet.h>  // Pour htonl/htons (Linux/macOS)

#include <cstring>
#include <string>

#include "../gauge/GaugeField.h"
#include "../mpi/MpiTopology.h"
// Inversion d'octets pour les doubles (Standard ILDG)
inline void swap_endian_64(void* ptr) {
    uint64_t* val = reinterpret_cast<uint64_t*>(ptr);
    *val = __builtin_bswap64(*val);
}

// Structure fixe de 144 octets pour LIME
struct LimeHeader {
    uint32_t magic;    // 0x454d494c (LIME en ASCII)
    uint16_t version;  // 1
    uint16_t mb_me;    // Message Begin / Message End flags
    uint64_t data_length;
    char type[128];  // Nom du record (ex: "ildg-binary-data")
};

std::string generate_ildg_xml(int L_glob);
void write_lime_header(MPI_File& fh, MPI_Offset offset, const std::string& type, uint64_t len,
                       bool mb, bool me);
void save_ildg_clime(const std::string& filename, const GaugeField& field, const GeometryCB& geo,
                     const mpi::MpiTopology& topo);
void read_ildg_clime(const std::string& filename, GaugeField& field, const GeometryCB& geo,
                     mpi::MpiTopology& topo);
