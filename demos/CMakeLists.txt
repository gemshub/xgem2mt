
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/gem2mt-data" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/nodearray-data" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/demo_nodearray_files.py" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/demo_nodearray_strings.py" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/demo_gem2mt_string.py" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/demo_gem2mt.py" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

set(DEMO_FILES
        demo_gem2mt
        demo_gem2mt_string
        demo_nodearray_files
        demo_nodearray_strings  )

foreach (exampleName ${DEMO_FILES})
    add_executable(${exampleName} ${exampleName}.cpp)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
          target_link_libraries(${exampleName}
                PUBLIC ${PROJECT_NAME}-static
                PUBLIC pthread )
   elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
          target_link_libraries( ${exampleName}
                PUBLIC ${PROJECT_NAME}
                PUBLIC pthread )
    else()
          target_link_libraries( ${exampleName}
               PUBLIC ${PROJECT_NAME}
               PUBLIC pthread )
    endif()
endforeach()
