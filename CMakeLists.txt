cmake_minimum_required(VERSION 3.10)
project(ECMC_QCB)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED) # Décommente si besoin
include(FetchContent)
#Lime
#FetchContent_Declare(
#        lime
#        GIT_REPOSITORY https://github.com/usqcd-software/c-lime.git
#        GIT_TAG        master  # Ou un tag spécifique comme v1.3.2
#)
#FetchContent_MakeAvailable(lime)
#Eigen
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found. Download via FetchContent...")

    FetchContent_Declare(
            eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(eigen)
else()
    message(STATUS "Eigen3 found on the system.")
endif()

#===ECMC CB===
add_executable(gauge_ecmc_cb src/bin/gauge_ecmc_cb.cpp
        src/gauge/GaugeField.cpp
        src/su3/utils.cpp
        src/geometry/GeometryCB.cpp
        src/geometry/types.h
        src/mpi/MpiTopology.cpp
        src/mpi/HalosExchange.cpp
)


# Configuration des options selon le compilateur détecté
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # Options spécifiques à icpx (OneAPI)
    target_compile_options(gauge_ecmc_cb PRIVATE
            -O2
            #-xHost            # Remplace -march=native pour Intel
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(gauge_ecmc_cb PRIVATE
            -O2
            #-ffast-math
            #-march=native
            -Wall -Wextra -Wpedantic
    )
endif()

target_link_libraries(gauge_ecmc_cb PRIVATE MPI::MPI_CXX Eigen3::Eigen)

