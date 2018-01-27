# CLion IDE project
#set(ENV{Path} "")
cmake_minimum_required(VERSION 3.5)
project(cec)

set(LINUX_R_INCLUDE_PATH "/usr/share/R/include")
set(LINUX_R "R")
set(LINUX_CEC_LIB "/home/konrad/R/x86_64-pc-linux-gnu-library/3.4/CEC/libs/CEC.so")

set(WINDOWS_R_INCLUDE_PATH "K:/R/R-3.4.3/include")
set(WINDOWS_R_LINK_DIRECTORY_PATH "K:/R/R-3.4.3/bin/x64")
set(WINDOWS_R "K:/R/R-3.4.3/bin/x64/R.exe")
set(WINDOWS_CEC_LIB "K:/R/R-3.4.3/library/CEC/libs/x64/CEC.dll")

if (WIN32)
    include_directories(${WINDOWS_R_INCLUDE_PATH})
    link_directories(${WINDOWS_R_LINK_DIRECTORY_PATH})
    link_libraries(m R Rlapack Rblas)
    set(R ${WINDOWS_R})
    set(CEC_LIB ${WINDOWS_CEC_LIB})
else (WIN32)
    include_directories(${LINUX_R_INCLUDE_PATH})
    link_libraries(m R lapack blas)
    set(R ${LINUX_R})
    set(CEC_LIB ${LINUX_CEC_LIB})
endif (WIN32)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(SOURCE_FILES
        src/cec_r.cpp
        src/cec_r.h
        src/cec_starter.cpp
        src/cec_starter.h
        src/cluster.h
        src/common.h
        src/exceptions.h
        src/init.cpp
        src/init.h
        src/m.h
        src/models/all.h
        src/models/cov_utils.cpp
        src/models/cov_utils.h
        src/models/covariance.h
        src/models/diagonal.h
        src/models/eigenvalues.h
        src/models/fixed_radius.h
        src/models/model.h
        src/models/spherical.h
        src/parallel_starter.cpp
        src/parallel_starter.h
        src/params.cpp
        src/params.h
        src/r_ext_ptr.h
        src/r_params.cpp
        src/r_params.h
        src/r_result.cpp
        src/r_result.h
        src/r_utils.h
        src/random.cpp
        src/random.h
        src/split_starter.cpp
        src/split_starter.h
        src/starter.cpp
        src/starter.h
        src/variable_starter.cpp
        src/variable_starter.h
        src/vec.h)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=core2")

set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -pipe -std=c++11 -pedantic -Wall")

set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg")

add_library(cec SHARED ${SOURCE_FILES})
#add_executable(cec_test ${SOURCE_FILES} src/test.cpp)

set(MAKEVARS_PKG_LIBS "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)")

file(WRITE src/Makevars "CXX_STD = CXX11" "\n")
file(APPEND src/Makevars "PKG_LIBS = " ${MAKEVARS_PKG_LIBS} "\n")

string(REPLACE ";" " " SOURCE_FILES_STRING "${SOURCE_FILES}")
string(REPLACE "src/" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REGEX REPLACE "[^ ]+\\.h" "" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REPLACE ".cpp" ".o" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REPLACE ".c" ".o" SOURCE_FILES_STRING ${SOURCE_FILES_STRING})
string(REGEX REPLACE " +" " " SOURCE_FILES_STRING ${SOURCE_FILES_STRING})

file(APPEND src/Makevars "OBJECTS = " ${SOURCE_FILES_STRING} "\n")

add_custom_target(install_CEC
        COMMAND ${R} CMD INSTALL --clean ${CMAKE_SOURCE_DIR}
        )

add_custom_target(dist_lib
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:cec> ${CEC_LIB}
        )

add_dependencies(dist_lib cec)
