# Collect all hpp files from the current directory
file(GLOB HPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Collect all cpp files from the current directory
file(GLOB CPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# The name of the shared and static libraries
set(XGEM2MT_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(XGEM2MT_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Build a xGEM2MT C++ library
add_library(xGEM2MT  ${HPP_FILES} ${CPP_FILES})

# Add an alias xGEM2MT::xGEM2MT to the target library xGEM2MT
add_library(xGEM2MT::xGEM2MT ALIAS xGEM2MT)

# Add the include paths to the target
target_include_directories(xGEM2MT
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Create install target for library
install(TARGETS xGEM2MT
    EXPORT xGEM2MTTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries)

# Install debug symbols
if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:xGEM2MT>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL
    )
endif()

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/xGEM2MT
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")

