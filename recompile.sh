#!/bin/bash
source modules_load.sh
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx
make -j $(nproc)
cd ..