//
// Created by ozdalkiran-l on 1/14/26.
//

#include <fstream>
#include <iostream>
#include <filesystem>
#include <map>
#include <string>
#include <sstream>
extern "C" {
    #include <lime.h>
}
#include "io.h"

namespace fs = std::filesystem;

//Saves a vector of doubles in ../data/filename.txt
void io::save_double(const std::vector<double> &data, const std::string &filename, int precision) {
    //Create a data folder if doesn't exists
    fs::path dir("data");

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir/ (filename+".txt");

    int counter = 1;
    while (fs::exists(filepath)) {
        //new name generation if file already exists
        std::string new_name = filename + "_" + std::to_string(counter) + ".txt";
        filepath = dir / new_name;
        counter++;
    }



    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::setprecision(precision);
    for (const double &x : data) {
        file << x << "\n";
    }
    file.close();
    std::cout << "\nData written in " << filepath << "\n";
}

//Saves a vector of doubles in ../data/filename and the run parameters in ../data/filename_params
void io::save_double_params(const std::vector<double> &data, const RunParams &params, const std::string &filename,
    int precision) {
    io::save_double(data, filename, precision);
    //No need to create data folder
    fs::path dir("data");

    //We create the params file
    std::string param_filename = filename + "_params.txt";
    fs::path filepath = dir/param_filename;

    int counter = 1;
    while (fs::exists(filepath)) {
        //new name generation if file already exists
        std::string new_name = filename + "_" + std::to_string(counter) + "_params.txt";
        filepath = dir / new_name;
        counter++;
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    //We write the run params
    file << "L_tot: " << params.L_core*params.n_core_dims << "\n"
            << "L_core: " << params.L_core << "\n"
            << "n_core_dims: " << params.n_core_dims << "\n"
            << "L_shift: " << params.L_shift << "\n"
            << "n_shift: " << params.n_shift << "\n"
            << "Seed: " << params.seed << "\n"
            << "beta: " << params.ecmc_params.beta << "\n"
            << "N_samples: " << params.ecmc_params.N_samples << "\n"
            << "param_theta_sample: " << params.ecmc_params.param_theta_sample << "\n"
            << "param_theta_refresh: " << params.ecmc_params.param_theta_refresh << "\n"
            << "poisson: " << (params.ecmc_params.poisson ? "true" : "false") << "\n"
            << "epsilon_set: " << params.ecmc_params.epsilon_set;
    file.close();
}

//Utilitary function to trim the spaces
std::string io::trim(const std::string& s) {
    size_t first = s.find_first_not_of(" \t");
    if (first == std::string::npos) return "";
    size_t last = s.find_last_not_of(" \t");
    return s.substr(first, (last - first + 1));
}

//Loads the parameters contained in filename into RunParams rp
void io::load_params(const std::string& filename, RunParams& rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Impossible d'ouvrir " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L_core"))      rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("L_shift"))     rp.L_shift = std::stoi(config["L_shift"]);
    if (config.count("n_shift"))     rp.n_shift = std::stoi(config["n_shift"]);
    if (config.count("stype_pos"))   rp.stype_pos = (config["stype_pos"] == "true");
    if (config.count("seed"))      rp.seed= std::stoi(config["seed"]);

    //ECMC params
    if (config.count("beta"))                rp.ecmc_params.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.ecmc_params.N_samples = std::stoi(config["N_samples"]);
    if (config.count("param_theta_sample"))  rp.ecmc_params.param_theta_sample = std::stod(config["param_theta_sample"]);
    if (config.count("param_theta_refresh")) rp.ecmc_params.param_theta_refresh = std::stod(config["param_theta_refresh"]);
    if (config.count("poisson"))             rp.ecmc_params.poisson = (config["poisson"] == "true");
    if (config.count("epsilon_set"))         rp.ecmc_params.epsilon_set = std::stod(config["epsilon_set"]);
}

//Loads the parameters contained in filename into RunParams rp
void io::load_params(const std::string& filename, RunParamsCB& rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Impossible d'ouvrir " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L_core"))      rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("n_shift"))     rp.n_shift = std::stoi(config["n_shift"]);
    if (config.count("stype_pos"))   rp.stype_pos = (config["stype_pos"] == "true");
    if (config.count("seed"))      rp.seed= std::stoi(config["seed"]);

    //ECMC params
    if (config.count("beta"))                rp.ecmc_params.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.ecmc_params.N_samples = std::stoi(config["N_samples"]);
    if (config.count("param_theta_sample"))  rp.ecmc_params.param_theta_sample = std::stod(config["param_theta_sample"]);
    if (config.count("param_theta_refresh")) rp.ecmc_params.param_theta_refresh = std::stod(config["param_theta_refresh"]);
    if (config.count("poisson"))             rp.ecmc_params.poisson = (config["poisson"] == "true");
    if (config.count("epsilon_set"))         rp.ecmc_params.epsilon_set = std::stod(config["epsilon_set"]);
}
//Loads the parameters for single core generation in filename into RunParamsSC rp
void io::load_params(const std::string &filename, RunParamsSC &rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open file " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L"))           rp.L = std::stoi(config["L"]);
    if (config.count("T"))           rp.T = std::stoi(config["T"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed"))  rp.seed = std::stoi(config["seed"]);


    //ECMC params
    if (config.count("beta"))                rp.ecmc_params.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.ecmc_params.N_samples = std::stoi(config["N_samples"]);
    if (config.count("param_theta_sample"))  rp.ecmc_params.param_theta_sample = std::stod(config["param_theta_sample"]);
    if (config.count("param_theta_refresh")) rp.ecmc_params.param_theta_refresh = std::stod(config["param_theta_refresh"]);
    if (config.count("poisson"))             rp.ecmc_params.poisson = (config["poisson"] == "true");
    if (config.count("epsilon_set"))         rp.ecmc_params.epsilon_set = std::stod(config["epsilon_set"]);
}

//Loads the parameters for Metropolis generations from input file
void io::load_params(const std::string &filename, RunParamsMetro &rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open file " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L"))           rp.L = std::stoi(config["L"]);
    if (config.count("T"))           rp.T = std::stoi(config["T"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed"))  rp.seed = std::stoi(config["seed"]);

    //Metropolis params
    if (config.count("beta"))                rp.mp.beta = std::stod(config["beta"]);
    if (config.count("epsilon"))         rp.mp.epsilon= std::stod(config["epsilon"]);
    if (config.count("N_samples"))           rp.mp.N_samples = std::stoi(config["N_samples"]);
    if (config.count("N_sweeps_meas"))           rp.mp.N_sweeps_meas = std::stoi(config["N_sweeps_meas"]);
    if (config.count("N_hits"))           rp.mp.N_hits = std::stoi(config["N_hits"]);
    if (config.count("N_burnin"))           rp.mp.N_burnin= std::stoi(config["N_burnin"]);
    if (config.count("N_set"))           rp.mp.N_set= std::stoi(config["N_set"]);
}

void io::load_params(const std::string &filename, RunParamsHb &rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open file " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L"))           rp.L = std::stoi(config["L"]);
    if (config.count("T"))           rp.T = std::stoi(config["T"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed"))  rp.seed = std::stoi(config["seed"]);

    //Hb params
    if (config.count("beta"))                rp.hp.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.hp.N_samples = std::stoi(config["N_samples"]);
    if (config.count("N_sweeps"))           rp.hp.N_sweeps = std::stoi(config["N_sweeps"]);
    if (config.count("N_hits"))           rp.hp.N_hits = std::stoi(config["N_hits"]);
}

void io::load_params(const std::string &filename, RunParamsHbMPI &rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open file " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L_core"))      rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("L_shift"))     rp.L_shift = std::stoi(config["L_shift"]);
    if (config.count("n_shift"))     rp.n_shift = std::stoi(config["n_shift"]);
    if (config.count("stype_pos"))   rp.stype_pos = (config["stype_pos"] == "true");
    if (config.count("seed"))      rp.seed= std::stoi(config["seed"]);

    //Hb params
    if (config.count("beta"))                rp.hp.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.hp.N_samples = std::stoi(config["N_samples"]);
    if (config.count("N_sweeps"))           rp.hp.N_sweeps = std::stoi(config["N_sweeps"]);
    if (config.count("N_hits"))           rp.hp.N_hits = std::stoi(config["N_hits"]);
}

void io::load_params(const std::string &filename, RunParamsECB& rp){
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Impossible d'ouvrir " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L_core"))      rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed"))      rp.seed= std::stoi(config["seed"]);
    if (config.count("N_switch_eo"))      rp.N_switch_eo= std::stoi(config["N_switch_eo"]);
    if (config.count("N_shift"))      rp.N_shift= std::stoi(config["N_shift"]);

    //ECMC params
    if (config.count("beta"))                rp.ecmc_params.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.ecmc_params.N_samples = std::stoi(config["N_samples"]);
    if (config.count("param_theta_sample"))  rp.ecmc_params.param_theta_sample = std::stod(config["param_theta_sample"]);
    if (config.count("param_theta_refresh")) rp.ecmc_params.param_theta_refresh = std::stod(config["param_theta_refresh"]);
    if (config.count("poisson"))             rp.ecmc_params.poisson = (config["poisson"] == "true");
    if (config.count("epsilon_set"))         rp.ecmc_params.epsilon_set = std::stod(config["epsilon_set"]);
}

void io::load_params(const std::string &filename, RunParamsHbCB& rp){
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open file " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        //Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    //Lattice params
    if (config.count("L_core"))      rp.L_core = std::stoi(config["L_core"]);
    if (config.count("n_core_dims")) rp.n_core_dims = std::stoi(config["n_core_dims"]);
    if (config.count("cold_start"))  rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed"))      rp.seed= std::stoi(config["seed"]);
    if (config.count("N_switch_eo"))      rp.N_switch_eo= std::stoi(config["N_switch_eo"]);
    if (config.count("N_shift"))      rp.N_shift= std::stoi(config["N_shift"]);

    //Hb params
    if (config.count("beta"))                rp.hp.beta = std::stod(config["beta"]);
    if (config.count("N_samples"))           rp.hp.N_samples = std::stoi(config["N_samples"]);
    if (config.count("N_sweeps"))           rp.hp.N_sweeps = std::stoi(config["N_sweeps"]);
    if (config.count("N_hits"))           rp.hp.N_hits = std::stoi(config["N_hits"]);
}

std::string io::ildg::generate_ildg_xml(int lx, int ly, int lz, int lt, int precision) {
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\" "
        << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        << "xsi:schemaLocation=\"http://www.lqcd.org/ildg http://www.lqcd.org/ildg/ildg-format.xsd\">\n"
        << "  <version>1.1</version>\n"
        << "  <field>su3gauge</field>\n"
        << "  <precision>" << precision << "</precision>\n"
        << "  <lx>" << lx << "</lx>\n"
        << "  <ly>" << ly << "</ly>\n"
        << "  <lz>" << lz << "</lz>\n"
        << "  <lt>" << lt << "</lt>\n"
        << "</ildgFormat>";
    return oss.str();
}

void io::ildg::save_ildg(const GaugeField &field, const GeometryCB &geo, const std::string &filename) {
    //Filepath for conf
    fs::path dir("data");
    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir/ (filename+".ildg");
    //If another file with same name exists
    int counter = 1;
    while (fs::exists(filepath)) {
        //new name generation if file already exists
        std::string new_name = filename + "_" + std::to_string(counter) + ".ildg";
        filepath = dir / new_name;
        counter++;
    }

    // 1. Prepare XML (double precision)
    std::string xml_header = generate_ildg_xml(geo.L, geo.L, geo.L, geo.L); // Lx, Ly, Lz, Lt

    // 2. Open C-LIME writer
    FILE *fp = fopen(filepath.c_str(), "wb");
    if (!fp) throw std::runtime_error("Impossible de créer le fichier");
    LimeWriter *writer = limeCreateWriter(fp);

    // 3. Write the XML Record
    int MB_flag = 1, ME_flag = 0;
    n_uint64_t xml_len = xml_header.size();
    LimeRecordHeader *h_xml = limeCreateHeader(MB_flag, ME_flag, (char*)"ildg-format", xml_len);
    limeWriteRecordHeader(h_xml, writer);
    limeWriteRecordData((void*)xml_header.c_str(), &xml_len, writer);
    limeDestroyHeader(h_xml);

    // 4. Prepare and write the binary data
    // Physical volume V (no halos)
    n_uint64_t nbytes = (n_uint64_t) geo.V * 4 * 9 * sizeof(std::complex<double>);
    ME_flag = 1;
    LimeRecordHeader *h_bin = limeCreateHeader(0, ME_flag, (char*)"ildg-binary-data", nbytes);
    limeWriteRecordHeader(h_bin, writer);

    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            // Copy of the matrix
            SU3 matrix = field.view_link_const(site, mu);

            // Endianess inversion
            double* raw_data = reinterpret_cast<double*>(matrix.data());
            for (int i = 0; i < 18; i++) {
                uint64_t *u = reinterpret_cast<uint64_t*>(&raw_data[i]);
                *u = __builtin_bswap64(*u);
            }

            // Matrix writing (18 doubles)
            n_uint64_t mat_size = 18 * sizeof(double);
            limeWriteRecordData(raw_data, &mat_size, writer);
        }
    }

    limeDestroyHeader(h_bin);
    limeDestroyWriter(writer);
    fclose(fp);
    std::cout << "Configuration written in " << filepath << "\n";
}

void io::ildg::read_ildg(GaugeField &field, const GeometryCB &geo, const std::string &filename) {
    fs::path dir("data");
    fs::path filepath = dir/ (filename+".ildg");
    FILE *fp = fopen(filepath.c_str(), "rb");
    if (!fp) throw std::runtime_error("Impossible d'ouvrir le fichier : " + filename);

    std::cout << "Reading configuration from " << filepath << "\n";
    LimeReader *reader = limeCreateReader(fp);
    bool data_found = false;

    // Parcourir les records du fichier LIME
    while (limeReaderNextRecord(reader) == LIME_SUCCESS) {
        char* type = limeReaderType(reader);

        // On cherche uniquement le record binaire
        if (std::string(type) == "ildg-binary-data") {
            data_found = true;

            n_uint64_t nbytes = limeReaderBytes(reader);
            n_uint64_t expected_bytes = (n_uint64_t) geo.V * 4 * 9 * sizeof(std::complex<double>);

            if (nbytes != expected_bytes) {
                throw std::runtime_error("Taille du fichier ILDG incohérente avec le volume du GaugeField !");
            }

            // Lecture site par site, mu par mu
            for (size_t site = 0; site < geo.V; ++site) {
                for (int mu = 0; mu < 4; ++mu) {

                    // On lit 18 doubles (9 complexes)
                    double buffer[18];
                    n_uint64_t to_read = 18 * sizeof(double);
                    limeReaderReadData(buffer, &to_read, reader);

                    // Conversion Big-Endian -> Little-Endian
                    for (int i = 0; i < 18; ++i) {
                        uint64_t *u = reinterpret_cast<uint64_t*>(&buffer[i]);
                        *u = __builtin_bswap64(*u);
                    }

                    // On "map" le buffer vers une matrice Eigen temporaire (Row-Major)
                    // Puis on l'affecte à notre view_link

                    // Méthode plus sûre pour copier les données vers votre vecteur 'links' via le Map
                    auto link = field.view_link(site, mu);
                    for(int i=0; i<3; ++i) {
                        for(int j=0; j<3; ++j) {
                            // ILDG stocke : Re(0,0), Im(0,0), Re(0,1), Im(0,1)...
                            double re = buffer[2 * (i * 3 + j)];
                            double im = buffer[2 * (i * 3 + j) + 1];
                            link(i, j) = std::complex<double>(re, im);
                        }
                    }
                }
            }
            break; // Données lues, on peut sortir
        }
    }

    if (!data_found) {
        throw std::runtime_error("Record 'ildg-binary-data' has not been found in the file.");
    }

    limeDestroyReader(reader);
    fclose(fp);
}
