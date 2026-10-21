# Copy and configure the doxygen configuration files to binary dir
configure_file(doxyfile-cpp      ${CMAKE_CURRENT_BINARY_DIR}/doxyfile-cpp @ONLY)

add_custom_target(docs ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/html/cpp"
    COMMAND doxygen doxyfile-cpp
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating API documentation with Doxygen..."
)
# Add target "docs-upload" for uploading the built documentation to reaktoro.org
#add_custom_target(docs-upload
#  COMMAND rsync -avz --delete ${CMAKE_CURRENT_BINARY_DIR}/html/ reaktoro:public_html
#  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

