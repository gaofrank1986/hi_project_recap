project(testf Fortran)
cmake_minimum_required(VERSION 3.4)
#set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# this is the directory, from which cmake was started, i.e. the top level source directory 
MESSAGE( STATUS "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR} )

# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

# the complete system name, e.g. "Linux-2.4.22", "FreeBSD-5.4-RELEASE" or "Windows 5.1" 
MESSAGE( STATUS "CMAKE_SYSTEM: " ${CMAKE_SYSTEM} )

Message(STATUS "Initial Fortran flags:" ${CMAKE_Fortran_FLAGS})
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  # something
  #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -warn all -ftrapuv -check all -debug all")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}   -traceback -warn all -ftrapuv")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}   -check all,noarg_temp_created")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -g -debug all")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  # something else
endif()

Message(STATUS "USER SET FORTRAN FLAGS:" ${CMAKE_Fortran_FLAGS})


ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(src/data)
ADD_SUBDIRECTORY(src/utils_motion)
ADD_SUBDIRECTORY(src/hi_intgl)
ADD_SUBDIRECTORY(src/lib_funcs)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

#add_library(hello SHARED ${SRC}/te.f90)
#add_executable(testf test.f90)

#target_link_libraries(testf hello)
#set_target_properties(hello PROPERTIES OUTPUT_NAME "hello")

