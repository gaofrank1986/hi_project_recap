
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)


add_library(data SHARED mesh.f90 wave.f90 motion.f90 time_mod.f90 potential.f90 data_wrapper.f90 misc_var_mod.f90 matrix_mod.f90)

#target_link_libraries(testf data)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

