# The path where cmake config files are installed
set(XGEM2MT_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ChemicalFun)

install(EXPORT xGEM2MTTargets
    FILE xGEM2MTTargets.cmake
    NAMESPACE xGEM2MT::
    DESTINATION ${XGEM2MT_INSTALL_CONFIGDIR}
    COMPONENT cmake)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/xGEM2MTConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/xGEM2MTConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/xGEM2MTConfig.cmake
    INSTALL_DESTINATION ${XGEM2MT_INSTALL_CONFIGDIR}
    PATH_VARS XGEM2MT_INSTALL_CONFIGDIR)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/xGEM2MTConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/xGEM2MTConfigVersion.cmake
    DESTINATION ${XGEM2MT_INSTALL_CONFIGDIR})
