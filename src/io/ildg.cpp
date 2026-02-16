#include "ildg.h"
#include <lime.h>

#include "../mpi/HalosExchange.h"

std::string generate_ildg_xml(int L_glob) {
    char buf[512];
    std::sprintf(buf,
                 "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                 "<ildgFormat><version>1.0</version><field>su3gauge</field>"
                 "<precision>64</precision><lx>%d</lx><ly>%d</ly><lz>%d</lz><lt>%d</lt>"
                 "</ildgFormat>",
                 L_glob, L_glob, L_glob, L_glob);
    return std::string(buf);
}

// Conversion en Big-Endian (requis par la norme ILDG)
inline double swap_double_endian(double val) {
    uint64_t bits;
    std::memcpy(&bits, &val, sizeof(bits));
    bits = __builtin_bswap64(bits);
    double swapped;
    std::memcpy(&swapped, &bits, sizeof(swapped));
    return swapped;
}

void save_ildg_clime(const std::string& filename, 
                     const GaugeField& field, 
                     const GeometryCB& geo, 
                     const mpi::MpiTopology& topo) 
{
    int L_global = geo.L_int * topo.n_core_dim;
    long nbytes_binary = (long)L_global * L_global * L_global * L_global * 72 * sizeof(double);
    size_t header_offset = 0;

    // --- ÉTAPE 1 : Le Rang 0 écrit les Headers LIME ---
    if (topo.rank == 0) {
        FILE *fp = fopen(filename.c_str(), "w");
        LimeWriter *w = limeCreateWriter(fp);

        // 1. Record XML (ildg-format)
        std::string xml = generate_ildg_xml(L_global);
        n_uint64_t xml_len = xml.size();
        LimeRecordHeader *h_xml = limeCreateHeader(1, 1, (char*)"ildg-format", xml.size());
        limeWriteRecordHeader(h_xml, w);
        limeWriteRecordData((void*)xml.c_str(), &xml_len, w);
        limeDestroyHeader(h_xml);

        // 2. Record Binaire (ildg-binary-data)
        // On écrit juste le header pour réserver l'espace des données
        LimeRecordHeader *h_bin = limeCreateHeader(0, 1, (char*)"ildg-binary-data", nbytes_binary);
        limeWriteRecordHeader(h_bin, w);
        
        // On récupère la position actuelle pour dire aux autres où commencer à écrire
        header_offset = ftell(fp);

        limeDestroyWriter(w);
        fclose(fp);
    }

    // --- ÉTAPE 2 : Synchronisation de l'offset ---
    MPI_Bcast(&header_offset, 1, MPI_OFFSET, 0, topo.cart_comm);

    // --- ÉTAPE 3 : Écriture MPI-I/O (Identique à la fonction précédente) ---
    // On réutilise la logique de local_buffer et de subarray
    
    std::vector<double> local_buffer(geo.V_int * 72);
    size_t buf_idx = 0;
    for (int t = 1; t <= geo.L_int; ++t) {
        for (int z = 1; z <= geo.L_int; ++z) {
            for (int y = 1; y <= geo.L_int; ++y) {
                for (int x = 1; x <= geo.L_int; ++x) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; ++mu) {
                        auto link = field.view_link_const(site, mu);
                        for (int i = 0; i < 3; ++i) {
                            for (int j = 0; j < 3; ++j) {
                                local_buffer[buf_idx++] = swap_double_endian(link(i, j).real());
                                local_buffer[buf_idx++] = swap_double_endian(link(i, j).imag());
                            }
                        }
                    }
                }
            }
        }
    }

    MPI_File fh;
    MPI_File_open(topo.cart_comm, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // Définition du type de donnée et du subarray (comme avant)
    MPI_Datatype site_type, file_type;
    MPI_Type_contiguous(72, MPI_DOUBLE, &site_type);
    MPI_Type_commit(&site_type);

    int g_sz[4] = {L_global, L_global, L_global, L_global};
    int l_sz[4] = {geo.L_int, geo.L_int, geo.L_int, geo.L_int};
    int start[4] = {topo.local_coords[3]*geo.L_int, topo.local_coords[2]*geo.L_int, 
                    topo.local_coords[1]*geo.L_int, topo.local_coords[0]*geo.L_int};

    MPI_Type_create_subarray(4, g_sz, l_sz, start, MPI_ORDER_C, site_type, &file_type);
    MPI_Type_commit(&file_type);

    // APPLICATION DE L'OFFSET LIME
    MPI_File_set_view(fh, (MPI_Offset)header_offset, site_type, file_type, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, local_buffer.data(), (int)geo.V_int, site_type, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
    MPI_Type_free(&file_type);
    MPI_Type_free(&site_type);
}

