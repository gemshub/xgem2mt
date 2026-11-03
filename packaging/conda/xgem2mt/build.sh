#!/bin/bash
set -eu
cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_PREFIX_PATH="$PREFIX" \
      -DCMAKE_CXX_FLAGS="${CXXFLAGS:-} -I$PREFIX/include/eigen3" \
      -DPYTHON_EXECUTABLE:FILEPATH="$PYTHON" \
      -DPython3_EXECUTABLE:FILEPATH="$PYTHON" \
      -DPython3_ROOT_DIR="$PREFIX" \
      -DPython3_FIND_FRAMEWORK=NEVER \
      -DPython3_FIND_STRATEGY=LOCATION \
      -DXGEM2MT_PYTHON_INSTALL_PREFIX="$PREFIX" \
      -DUSE_SPDLOG_PRECOMPILED=ON \
      -DXGEM2MT_BUILD_DEMOS=OFF \
      -DXGEM2MT_BUILD_TESTS=OFF
cmake --build build -j"${CPU_COUNT}"
cmake --install build
