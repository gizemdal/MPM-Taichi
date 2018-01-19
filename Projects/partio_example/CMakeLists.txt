#Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_cispba_executable(partio_example main.cpp)

target_link_libraries(partio_example partio)
