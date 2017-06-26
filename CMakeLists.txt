# CLion IDE project
cmake_minimum_required(VERSION 3.5)
project(cec)
set(CMAKE_C_STANDARD 99)

include_directories($ENV{R_INCLUDE_PATH})

set(SOURCE_FILES
        src/alloc.h
        src/array.c
        src/array.h
        src/cec.c
        src/cec.h
        src/cec_context.c
        src/cec_context.h
        src/cec_params.h
        src/cec_params_r.c
        src/cec_params_r.h
        src/cec_r.c
        src/cec_r.h
        src/cec_r_utils.c
        src/cec_r_utils.h
        src/cec_starter_omp.c
        src/cec_starter_omp.h
        src/centers_init.c
        src/centers_init.h
        src/centers_init_r.c
        src/centers_init_r.h
        src/cov_utils.c
        src/cov_utils.h
        src/cross_entropy_context.h
        src/error_r.c
        src/error_r.h
        src/errors.h
        src/init_rand.c
        src/init_rand.h
        src/init_utils_r.c
        src/init_utils_r.h
        src/kmeanspp.c
        src/kmeanspp.h
        src/matrix.c
        src/matrix.h
        src/mem_mg.c
        src/mem_mg.h
        src/model.c
        src/model.h
        src/noret.h
        src/rand.c
        src/rand.h
        src/result_r.c
        src/result_r.h
        src/models/all.c
        src/models/covariance.c
        src/models/diagonal.c
        src/models/eigenvalues.c
        src/models/energy.h
        src/models/fixedr.c
        src/models/models.c
        src/models/models.h
        src/models/spherical.c
        )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -DR_ALLOC -DR_RAND -Wall -Werror")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z defs")

link_libraries(m R lapack blas)

add_library(cec SHARED ${SOURCE_FILES})

set(MAKEVARS_PKG_CFLAGS "-DR_ALLOC -DR_RAND $(SHLIB_OPENMP_CFLAGS)")
set(MAKEVARS_PKG_LIBS "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CFLAGS)")

file(WRITE src/Makevars "PKG_CFLAGS = " ${MAKEVARS_PKG_CFLAGS} "\n")
file(APPEND src/Makevars "PKG_LIBS = " ${MAKEVARS_PKG_LIBS} "\n")

string (REPLACE ";" " " SOURCE_FILES_STRING "${SOURCE_FILES}")
string (REPLACE "src/" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string (REGEX REPLACE "[^ ]+\\.h" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string (REPLACE ".c" ".o" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string (REGEX REPLACE " +" " " SOURCE_FILES_STRING ${SOURCE_FILES_STRING})

file(APPEND src/Makevars "OBJECTS = " ${SOURCE_FILES_STRING} "\n")

add_custom_target(install_CEC
        COMMAND R CMD INSTALL --clean ${CMAKE_SOURCE_DIR}
        )

add_custom_target(dist_lib
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:cec> $ENV{R_CEC_LIB_PATH}
        )
