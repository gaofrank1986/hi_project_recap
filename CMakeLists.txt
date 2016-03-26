project(testf Fortran)

#set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(src/data)
ADD_SUBDIRECTORY(src/property)
ADD_SUBDIRECTORY(src/hi_intgl)
ADD_SUBDIRECTORY(src/lib_funcs)


#add_library(hello SHARED ${SRC}/te.f90)
#add_executable(testf test.f90)

#target_link_libraries(testf hello)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

