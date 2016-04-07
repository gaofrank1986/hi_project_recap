
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_library(utils SHARED gaussian.f90 free_term.f90 grad.f90)
#add_executable(testf test.f90)

ADD_SUBDIRECTORY(./hi_intgl)

target_link_libraries(utils func const hi_mod)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

