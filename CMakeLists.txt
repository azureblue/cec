# CLion IDE project
cmake_minimum_required(VERSION 3.5)
project(cec)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

include_directories(/usr/share/R/include)

set(SOURCE_FILES
        src/cec_params.cpp
        src/cec_params.h
        src/cec_params_r.cpp
        src/cec_params_r.h
        src/cec_starter.cpp
        src/cec_starter.h
        src/r_utils.cpp
        src/r_utils.h
        src/cluster.h
        src/exceptions.h
        src/single_start_input.h
        src/models/cov_utils.h
        src/models/cov_utils.cpp
        src/models/model.h
        src/models/all.h
        src/vec.h
        src/result_r.h
        src/result_r.cpp
        src/cec_r.cpp
        src/cec_r.h
        src/init.h)


set(CMAKE_C_FLAGS "-O0 -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -DR_ALLOC -DR_RAND -Wall")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z defs -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z defs -fopenmp")

link_libraries(m R lapack blas)

add_library(cec SHARED ${SOURCE_FILES} src/exceptions.h src/models/cov_utils.cpp)
#add_executable(cec_test ${SOURCE_FILES} src/vec_test.cpp)

set(MAKEVARS_PKG_CFLAGS "-DR_ALLOC -DR_RAND $(SHLIB_OPENMP_CFLAGS)")
set(MAKEVARS_PKG_LIBS "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CFLAGS)")

file(WRITE src/Makevars "PKG_CFLAGS = " ${MAKEVARS_PKG_CFLAGS} "\n")
file(APPEND src/Makevars "PKG_LIBS = " ${MAKEVARS_PKG_LIBS} "\n")

string (REPLACE ";" " " SOURCE_FILES_STRING "${SOURCE_FILES}")
string (REPLACE "src/" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string (REGEX REPLACE "[^ ]+\\.h" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string (REPLACE ".cpp" ".o" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string (REPLACE ".c" ".o" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string (REGEX REPLACE " +" " " SOURCE_FILES_STRING ${SOURCE_FILES_STRING})

file(APPEND src/Makevars "OBJECTS = " ${SOURCE_FILES_STRING} "\n")

add_custom_target(install_CEC
        COMMAND R CMD INSTALL --clean ${CMAKE_SOURCE_DIR}
        )

add_custom_target(dist_lib
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:cec> $ENV{R_CEC_LIB_PATH}
        )

add_dependencies(dist_lib cec)