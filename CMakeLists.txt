# CLion IDE project
set(ENV{Path} "" )
cmake_minimum_required(VERSION 3.5)
project(cec)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

include_directories(/usr/share/R/include)

set(SOURCE_FILES
        src/params.cpp
        src/params.h
        src/r_params.cpp
        src/r_params.h
        src/starter.cpp
        src/starter.h
        src/r_ext_ptr.h
        src/r_utils.h
        src/cluster.h
        src/exceptions.h
        src/models/cov_utils.h
        src/models/cov_utils.cpp
        src/models/model.h
        src/models/all.h
        src/models/diagonal.h
        src/models/spherical.h
        src/models/fixed_radius.h
        src/models/covariance.h
        src/models/eigenvalues.h
        src/vec.h
        src/r_result.h
        src/r_result.cpp
        src/cec_r.cpp
        src/cec_r.h
        src/init.h src/models/constants.h)


link_libraries(m R lapack blas)
set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -std=c++0x -pedantic -Wall")
add_library(cec SHARED ${SOURCE_FILES})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
add_executable(cec_test ${SOURCE_FILES} src/test.cpp src/test.cpp)

set(MAKEVARS_PKG_CFLAGS "-DR_ALLOC -DR_RAND $(SHLIB_OPENMP_CFLAGS)")
set(MAKEVARS_PKG_LIBS "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CFLAGS)")

file(WRITE src/Makevars "CXX_STD = CXX11" "\n")
file(APPEND src/Makevars "PKG_CFLAGS = " ${MAKEVARS_PKG_CFLAGS} "\n")
file(APPEND src/Makevars "PKG_LIBS = " ${MAKEVARS_PKG_LIBS} "\n")

string(REPLACE ";" " " SOURCE_FILES_STRING "${SOURCE_FILES}")
string(REPLACE "src/" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REGEX REPLACE "[^ ]+\\.h" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REPLACE ".cpp" ".o" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REPLACE ".c" ".o" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REGEX REPLACE " +" " " SOURCE_FILES_STRING ${SOURCE_FILES_STRING})

file(APPEND src/Makevars "OBJECTS = " ${SOURCE_FILES_STRING} "\n")

add_custom_target(install_CEC
        COMMAND R CMD INSTALL --clean ${CMAKE_SOURCE_DIR}
        )

add_custom_target(dist_lib
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:cec> /home/konrad/R/x86_64-pc-linux-gnu-library/3.4/CEC/libs/CEC.so
        )

add_dependencies(dist_lib cec)
