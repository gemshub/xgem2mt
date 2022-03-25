# Find pybind11
#find_package(pybind11 REQUIRED)

# Build the Python wrappers of xGEM2MT library
add_subdirectory(PyxGEM2MT)

# Build the Python package for xGEM2MT
add_subdirectory(xGEM2MT)
