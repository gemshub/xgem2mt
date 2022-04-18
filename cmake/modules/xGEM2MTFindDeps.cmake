# Find Eigen library and include its directory
#find_package(Eigen3 3.3.7 REQUIRED)
#message(STATUS "Found Eigen3 v${Eigen3_VERSION}: ${Eigen3_INCLUDE_DIRS}")


# Find pybind11 library (if needed)
if(XGEM2MT_BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11 - the Python module `xGEM2MT` will not be built.")
        set(XGEM2MT_BUILD_PYTHON OFF)
    else()
        message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    endif()
endif()
