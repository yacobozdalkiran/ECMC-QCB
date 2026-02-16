#include "ildg.h"

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

void write_lime_header(MPI_File& fh, MPI_Offset offset, const std::string& type, uint64_t len,
                       bool mb, bool me) {
    LimeHeader h;
    std::memset(&h, 0, 144);
    h.magic = htonl(0x454d494c);
    h.version = htons(1);
    h.mb_me = htons((mb ? 1 : 0) | (me ? 2 : 0));
    h.data_length = __builtin_bswap64(len);
    std::strncpy(h.type, type.c_str(), 127);

    MPI_File_write_at(fh, offset, &h, 144, MPI_BYTE, MPI_STATUS_IGNORE);
}

void save_configuration_lime(const std::string& filename, const GaugeField& field,
                             const GeometryCB& geo, const mpi::MpiTopology& topo) {
    int L_glob = geo.L_int * topo.n_core_dim;
    std::string xml = generate_ildg_xml(L_glob);
    uint64_t n_links_glob = (uint64_t)L_glob * L_glob * L_glob * L_glob * 4;
    uint64_t binary_size = n_links_glob * 18 * sizeof(double);

    // Offsets LIME : Header1(144) + XML + Header2(144) + Data
    MPI_Offset xml_off = 144;
    MPI_Offset header2_off = xml_off + xml.size();
    MPI_Offset data_off = header2_off + 144;

    MPI_File fh;
    MPI_File_open(topo.cart_comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

    if (topo.rank == 0) {
        write_lime_header(fh, 0, "ildg-format", xml.size(), true, false);
        MPI_File_write_at(fh, xml_off, xml.data(), xml.size(), MPI_BYTE, MPI_STATUS_IGNORE);
        write_lime_header(fh, header2_off, "ildg-binary-data", binary_size, false, true);
    }

    // Vue MPI : on définit le sous-cube local dans le cube global (T, Z, Y, X, matrices_18_doubles)
    int g_sizes[5] = {L_glob, L_glob, L_glob, L_glob, 18};
    int l_sizes[5] = {geo.L_int, geo.L_int, geo.L_int, geo.L_int, 18};
    int starts[5] = {topo.local_coords[3] * geo.L_int, topo.local_coords[2] * geo.L_int,
                     topo.local_coords[1] * geo.L_int, topo.local_coords[0] * geo.L_int, 0};

    MPI_Datatype file_view;
    MPI_Type_create_subarray(5, g_sizes, l_sizes, starts, MPI_ORDER_C, MPI_DOUBLE, &file_view);
    MPI_Type_commit(&file_view);

    // Préparation du buffer local (Extraction du coeur sans les halos + Swap Endian)
    std::vector<double> buffer;
    buffer.reserve(geo.V_int * 4 * 18);
    for (int t = 1; t <= geo.L_int; ++t) {
        for (int z = 1; z <= geo.L_int; ++z) {
            for (int y = 1; y <= geo.L_int; ++y) {
                for (int x = 1; x <= geo.L_int; ++x) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; ++mu) {
                        auto link = field.view_link_const(site, mu);
                        for (int i = 0; i < 3; ++i)
                            for (int j = 0; j < 3; ++j) {
                                double re = link(i, j).real(), im = link(i, j).imag();
                                swap_endian_64(&re);
                                swap_endian_64(&im);
                                buffer.push_back(re);
                                buffer.push_back(im);
                            }
                    }
                }
            }
        }
    }

    MPI_File_set_view(fh, data_off, MPI_DOUBLE, file_view, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, buffer.data(), buffer.size(), MPI_DOUBLE, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
    MPI_Type_free(&file_view);
}

void load_configuration_lime(const std::string& filename, GaugeField& field, const GeometryCB& geo,
                             mpi::MpiTopology& topo) {
    int L_glob = geo.L_int * topo.n_core_dim;
    std::string xml = generate_ildg_xml(L_glob);
    MPI_Offset data_off = 144 + xml.size() + 144;

    MPI_File fh;
    if (MPI_File_open(topo.cart_comm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) !=
        MPI_SUCCESS) {
        if (topo.rank == 0) std::cerr << "Erreur : Impossible d'ouvrir " << filename << std::endl;
        return;
    }

    int g_sizes[5] = {L_glob, L_glob, L_glob, L_glob, 18};
    int l_sizes[5] = {geo.L_int, geo.L_int, geo.L_int, geo.L_int, 18};
    int starts[5] = {topo.local_coords[3] * geo.L_int, topo.local_coords[2] * geo.L_int,
                     topo.local_coords[1] * geo.L_int, topo.local_coords[0] * geo.L_int, 0};

    MPI_Datatype file_view;
    MPI_Type_create_subarray(5, g_sizes, l_sizes, starts, MPI_ORDER_C, MPI_DOUBLE, &file_view);
    MPI_Type_commit(&file_view);

    std::vector<double> buffer(geo.V_int * 4 * 18);
    MPI_File_set_view(fh, data_off, MPI_DOUBLE, file_view, "native", MPI_INFO_NULL);
    MPI_File_read_all(fh, buffer.data(), buffer.size(), MPI_DOUBLE, MPI_STATUS_IGNORE);

    // Reconstruction des matrices SU(3) dans GaugeField
    size_t idx = 0;
    for (int t = 1; t <= geo.L_int; ++t) {
        for (int z = 1; z <= geo.L_int; ++z) {
            for (int y = 1; y <= geo.L_int; ++y) {
                for (int x = 1; x <= geo.L_int; ++x) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; ++mu) {
                        SU3 mat;
                        for (int i = 0; i < 3; ++i)
                            for (int j = 0; j < 3; ++j) {
                                double re = buffer[idx++], im = buffer[idx++];
                                swap_endian_64(&re);
                                swap_endian_64(&im);
                                mat(i, j) = std::complex<double>(re, im);
                            }
                        field.view_link(site, mu) = mat;
                    }
                }
            }
        }
    }

    MPI_File_close(&fh);
    MPI_Type_free(&file_view);

    mpi::exchange::exchange_halos_cascade(field, geo, topo);
}
