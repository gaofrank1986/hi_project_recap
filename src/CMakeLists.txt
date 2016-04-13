
set(SRC_LIST  elem_proc.f90 tassbt.f90 tforce.f90  main.f90 timestep.f90 tass0_time.f90)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

#ADD_SUBDIRECTORY(${SRC} bin)


#add_library(hello SHARED ${SRC}/te.f90)
add_executable(testf ${SRC_LIST})

target_link_libraries(testf data utils func const)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

