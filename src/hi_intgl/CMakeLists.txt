set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
add_library(hi_mod SHARED hi_funcs.f90 hi_integral.f90)

target_link_libraries(hi_mod data)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

