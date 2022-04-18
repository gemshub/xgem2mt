# Collect the source files of PyxGEM2MT
file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Include the xGEM2MT/python directory
include_directories(${PROJECT_SOURCE_DIR}/python)
include_directories(${CONDA_AWARE_PREFIX}/include/eigen3)
include_directories("/usr/local/include/eigen3")

# Create the PyxGEM2MT library
pybind11_add_module(PyxGEM2MT ${SOURCE_FILES} NO_EXTRAS) # NO_EXTRAS used to decrease linking time, but causes larger module file size

# Link PyxGEM2MT against xGEM2MT and other library
target_link_libraries(PyxGEM2MT LINK_PRIVATE xGEM2MT::xGEM2MT GEMS3K-static GEM2MT-static)
