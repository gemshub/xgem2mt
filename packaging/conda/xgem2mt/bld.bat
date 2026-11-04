set "PYTHON_FWD=%PYTHON:\=/%"

cmake -S . -B build -G Ninja ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
      -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
      -DPYTHON_EXECUTABLE:FILEPATH="%PYTHON_FWD%" ^
      -DPython3_EXECUTABLE:FILEPATH="%PYTHON_FWD%" ^
      -DPython3_INCLUDE_DIR="%PYTHON_INCLUDE%" ^
      -DPython3_LIBRARY="%PYTHON_LIB%" ^
      -DXGEM2MT_PYTHON_INSTALL_PREFIX="%PREFIX%" ^
      -DUSE_SPDLOG_PRECOMPILED=ON ^
      -DXGEM2MT_BUILD_DEMOS=OFF ^
      -DXGEM2MT_BUILD_TESTS=OFF ^
      -DCMAKE_CXX_FLAGS="/utf-8 /I%LIBRARY_PREFIX%\include\eigen3"
if errorlevel 1 exit 1
cmake --build build -j 1
if errorlevel 1 exit 1
cmake --install build
if errorlevel 1 exit 1
