# Collect all hpp files from the current directory
file(GLOB HPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Collect all cpp files from the current directory
file(GLOB CPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# The name of the shared and static libraries
#set(XGEM2MT_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
#set(XGEM2MT_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Compile the C++ source code into object files
add_library(XGEM2MT_OBJECT OBJECT ${HPP_FILES} ${CPP_FILES})

# Check if a shared library is to be built
if(BUILD_SHARED_LIBS)
    # Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    add_library(xGEM2MT SHARED $<TARGET_OBJECTS:XGEM2MT_OBJECT>)
    add_library(xGEM2MT::xGEM2MT ALIAS xGEM2MT)
    target_link_libraries(xGEM2MT PUBLIC GEMS3K::GEMS3K GEMS3K::GEM2MT)
    install( TARGETS xGEM2MT
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
if(BUILD_STATIC_LIBS)
   add_library(xGEM2MT-static STATIC $<TARGET_OBJECTS:XGEM2MT_OBJECT>)
   add_library(xGEM2MT::xGEM2MT-static ALIAS xGEM2MT-static)
   target_link_libraries(xGEM2MT-static PUBLIC GEMS3K::GEMS3K-static GEMS3K::GEM2MT-static)
   install( TARGETS xGEM2MT-static
             EXPORT xGEM2MTTargets DESTINATION "lib" COMPONENT libraries)
endif()

# Install the header files preserving the directory hierarchy
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp" )

