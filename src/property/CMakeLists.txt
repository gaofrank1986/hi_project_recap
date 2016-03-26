
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_library(property SHARED gaussian.f90 bd_prop.f90)
#add_executable(testf test.f90)

target_link_libraries(property func)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

