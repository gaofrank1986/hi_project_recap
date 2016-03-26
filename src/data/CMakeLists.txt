
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)


add_library(data SHARED mesh.f90 wave.f90)

#target_link_libraries(testf data)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

