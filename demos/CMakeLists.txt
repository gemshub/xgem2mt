
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/gem2mt-data" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/nodearray-data" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

set(DEMO_FILES
        demo_gem2mt
        demo_nodearray_files
        demo_nodearray_strings  )

foreach (exampleName ${DEMO_FILES})
    add_executable(${exampleName} ${exampleName}.cpp)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
          target_link_libraries(${exampleName}
                PUBLIC ${PROJECT_NAME}
                PUBLIC GEMS3K-static
                PUBLIC GEM2MT-static
                PUBLIC pthread )
   elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
          target_link_libraries( ${exampleName}
                PUBLIC ${PROJECT_NAME}
                PUBLIC GEMS3K
                PUBLIC GEM2MT
                PUBLIC pthread )
    else()
          target_link_libraries( ${exampleName}
               PUBLIC ${PROJECT_NAME}
               PUBLIC GEMS3K
               PUBLIC GEM2MT
               PUBLIC pthread )
    endif()
endforeach()
