# Collect all hpp files from the current directory
file(GLOB HPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Collect all cpp files from the current directory
file(GLOB CPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Compile the C++ source code into object files
add_library(XGEM2MT_OBJECT OBJECT ${HPP_FILES} ${CPP_FILES})

# Compile with the spdlog usage requirements (SPDLOG_COMPILED_LIB, SPDLOG_SHARED_LIB, ...).
# Without them MSVC emits its own copies of the spdlog::logger functions, which then clash
# with the ones imported from spdlog.dll when the Python module is linked (LNK2005).
if(TARGET spdlog::spdlog)
    target_link_libraries(XGEM2MT_OBJECT PUBLIC spdlog::spdlog)
endif()

# Check if a shared library is to be built
if(XGEM2MT_BUILD_SHARED_LIBS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    add_library(xGEM2MT SHARED $<TARGET_OBJECTS:XGEM2MT_OBJECT>)
    add_library(xGEM2MT::xGEM2MT ALIAS xGEM2MT)
    target_link_libraries(xGEM2MT PUBLIC GEMS3K::GEMS3K)
    install(TARGETS xGEM2MT
            EXPORT xGEM2MTTargets DESTINATION "lib" COMPONENT libraries)

    # Install debug symbols
    if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:xGEM2MT>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL)
    endif()
endif()

# Check if a static library is to be built
if(XGEM2MT_BUILD_STATIC_LIBS)
   add_library(xGEM2MT-static STATIC $<TARGET_OBJECTS:XGEM2MT_OBJECT>)
   add_library(xGEM2MT::xGEM2MT-static ALIAS xGEM2MT-static)
   target_link_libraries(xGEM2MT-static PUBLIC GEMS3K::GEMS3K-static)
   install(TARGETS xGEM2MT-static
             EXPORT xGEM2MTTargets DESTINATION "lib" COMPONENT libraries)
endif()

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp" )

