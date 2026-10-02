message(STATUS "Building python")

# Configure the pyproject.toml script with the provided CMake variables
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml.in ${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml)

add_subdirectory(src/xgem2mt)

# Install where the python package is installed to CMAKE_INSTALL_PREFIX if not given
if(NOT DEFINED XGEM2MT_PYTHON_INSTALL_PREFIX)

# Install the xgem2mt python package using setuptools
install(CODE
"
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyxGEM2MT>
                ${CMAKE_CURRENT_BINARY_DIR}/src/xgem2mt/$<TARGET_FILE_NAME:PyxGEM2MT>
      )

   execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install ${CMAKE_CURRENT_BINARY_DIR} -vv
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
"
)

else()

# If the path is already in Windows format (with backslashes), it can't be added directly
# to the string below, otherwise CMake will later complain about "Invalid escape sequence".
file(TO_CMAKE_PATH "${XGEM2MT_PYTHON_INSTALL_PREFIX}" XGEM2MT_PYTHON_INSTALL_PREFIX)
message(STATUS "XGEM2MT_PYTHON_INSTALL_PREFIX: ${XGEM2MT_PYTHON_INSTALL_PREFIX}")

# Install the xgem2mt python package using setuptools
install(CODE
"
   file(TO_NATIVE_PATH \"${XGEM2MT_PYTHON_INSTALL_PREFIX}\" XGEM2MT_PYTHON_INSTALL_PREFIX_NATIVE)

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyxGEM2MT>
                ${CMAKE_CURRENT_BINARY_DIR}/src/xgem2mt/$<TARGET_FILE_NAME:PyxGEM2MT>
        )

    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install --prefix \${XGEM2MT_PYTHON_INSTALL_PREFIX_NATIVE} ${CMAKE_CURRENT_BINARY_DIR} -vv
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
"
)

endif()

