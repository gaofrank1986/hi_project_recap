set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
add_library(func SHARED constant.f90 kinds.f90 shape_func.f90 wave_func.f90 green_func.f90 solidagM.f90
    linalg.f90 tripole.f90 WFuncs1.f90)

target_link_libraries(func data)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

