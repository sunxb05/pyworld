cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(PElib LANGUAGES Fortran VERSION 1.6)

include(GNUInstallDirs)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR} CACHE PATH "Set the library build directory")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR} CACHE PATH "Set the archive build directory")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Set the Fortran module build directory")

option(ENABLE_PDE "Enable PDE" OFF)
option(ENABLE_MPI "Enable MPI" OFF)
option(USE_MPIF "Use mpif.h instead of mpi.mod" OFF)
option(ENABLE_COVERAGE "Enable code coverage (requires GCC)" OFF)

set(PElib_SOURCES
    ${PROJECT_SOURCE_DIR}/src/pelib_mpi.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_precision.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_constants.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_options.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_utils.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_pde.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_operators.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_multipole_tools.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_interaction_tools.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_analysis_tools.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_cavity_generators.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_potential_derivatives.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_induced_moments.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/pelib_blas_interfaces.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/pelib_lapack_interfaces.F90
    ${PROJECT_SOURCE_DIR}/src/interfaces/pelib_integral_interfaces.F90
    ${PROJECT_SOURCE_DIR}/src/fmm/fmm.F90
    ${PROJECT_SOURCE_DIR}/src/fmm/tensors.F90
    ${PROJECT_SOURCE_DIR}/src/fmm/tensors_recursive.F90
    ${PROJECT_SOURCE_DIR}/src/fmm/tensors_damped_erf.F90
    ${PROJECT_SOURCE_DIR}/src/fmm/tensors_damped_thole.F90
    ${PROJECT_SOURCE_DIR}/src/fmm/tensors_damped_amoeba.F90
    ${PROJECT_SOURCE_DIR}/src/pelib_options.F90
    ${PROJECT_SOURCE_DIR}/src/pelib.F90
    )

add_library(PElib ${PElib_SOURCES})

set(REAL_PRECISION "REAL64" CACHE STRING "Choose precision of reals.")
set_property(CACHE REAL_PRECISION PROPERTY STRINGS "REAL32" "REAL64" "REAL128")
if(${REAL_PRECISION} MATCHES "REAL32")
    message(STATUS "Using 32-bit real precision.")
    target_compile_definitions(PElib PRIVATE REAL32)
elseif(${REAL_PRECISION} MATCHES "REAL64")
    message(STATUS "Using 64-bit real precision.")
    target_compile_definitions(PElib PRIVATE REAL64)
elseif(${REAL_PRECISION} MATCHES "REAL128")
    message(STATUS "Using 128-bit real precision.")
    target_compile_definitions(PElib PRIVATE REAL128)
else()
    message(FATAL_ERROR "Unsupported real precision: ${REAL_PRECISION}" )
endif()

set(INTEGER_PRECISION "INT32" CACHE STRING "Choose precision of integers.")
set_property(CACHE INTEGER_PRECISION PROPERTY STRINGS "INT32" "INT64")
if(${INTEGER_PRECISION} MATCHES "INT32")
    message(STATUS "Using 32-bit integer precision.")
    target_compile_definitions(PElib PRIVATE INT32)
elseif(${INTEGER_PRECISION} MATCHES "INT64")
    message(STATUS "Using 64-bit integer precision.")
    target_compile_definitions(PElib PRIVATE INT64)
else()
    message(FATAL_ERROR "Unsupported integer precision: ${INTEGER_PRECISION}")
endif()

set(CMAKE_Fortran_FORMAT "FREE")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

message(STATUS "Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message(STATUS "System processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(REAL_PRECISION MATCHES "REAL64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
    elseif(REAL_PRECISION MATCHES "REAL128")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r16")
    endif()
    if(INTEGER_PRECISION MATCHES "INT64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
    endif()
    if(ENABLE_COVERAGE)
        message(FATAL_ERROR "Code coverage compilation requires GCC compiler.")
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -no-prec-div -fp-model fast=2 -ip -xHost")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -traceback -check all -debug all -ftrapuv -warn all")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
    if(REAL_PRECISION MATCHES "REAL64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
    elseif(REAL_PRECISION MATCHES "REAL128")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-16")
    endif()
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0.0)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()
    if(INTEGER_PRECISION MATCHES "INT64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
    endif()
    if(ENABLE_COVERAGE)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --coverage")
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -Ofast -march=native -mtune=native")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Og -fcheck=all -fbacktrace -Wall -Wextra -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wsurprising -Wimplicit-interface -Wimplicit-procedure -Wrealloc-lhs -Wrealloc-lhs-all")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "8.2.0")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -std=f2008ts")
    else()
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -std=f2018")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    if(NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
        if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "aarch64")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=large")
        else()
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
        endif()
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -Mipa=fast")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    if(REAL_PRECISION MATCHES "REAL64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
    elseif(REAL_PRECISION MATCHES "REAL128")
        message(WARNING "${REAL_PRECISION} not available through ${CMAKE_Fortran_COMPILER_ID}")
    endif()
    if(INTEGER_PRECISION MATCHES "INT64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
    endif()
    if(NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
        if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "aarch64")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=large")
        else()
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
        endif()
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -Mipa=fast")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qzerosize -qextname -qfree -qlanglvl=extended -qinit=f90ptr")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -qstrict -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -eZ")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")
else()
    message(WARNING "PElib has not been explicitly configured for ${CMAKE_Fortran_COMPILER}.")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
else()
    message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS}")
endif()

if(ENABLE_MPI)
    target_compile_definitions(PElib PRIVATE VAR_MPI)
    if(NOT USE_MPIF)
        target_compile_definitions(PElib PRIVATE USE_MPI_MOD_F90)
    endif()
endif()

if(ENABLE_MPI AND NOT DEFINED HOST_PROGRAM)
    find_package(MPI REQUIRED Fortran)
    target_compile_options(PElib PRIVATE ${MPI_Fortran_COMPILE_OPTIONS})
    target_compile_definitions(PElib PRIVATE ${MPI_Fortran_COMPILE_DEFINITIONS})
    target_include_directories(PElib PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
endif()

if(ENABLE_PDE)
    target_compile_definitions(PElib PRIVATE VAR_PDE)
endif()

if(ENABLE_PDE AND NOT DEFINED HOST_PROGRAM)
    find_package(HDF5 REQUIRED Fortran)
    target_include_directories(PElib PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
endif()

if(DEFINED HOST_PROGRAM)
    if("${INTEGRAL_LIBRARY}" STREQUAL "GEN1INT")
        target_compile_definitions(PElib PRIVATE GEN1INT_INTEGRALS)
        message(STATUS "Using Gen1Int integral library.")
    else()
        message(FATAL_ERROR "No integral library enabled. Options: GEN1INT")
    endif()
else()
    target_compile_definitions(PElib PRIVATE DUMMY_INTEGRALS)
    message(STATUS "Using dummy integral interfaces.")
endif()

if("${HOST_PROGRAM}" STREQUAL "DIRAC")
    target_compile_definitions(PElib PRIVATE PRG_DIRAC)
    message(STATUS "Host program set to DIRAC.")
elseif("${HOST_PROGRAM}" STREQUAL "DALTON")
    message(STATUS "Host program set to Dalton.")
endif()

if(DEFINED PARENT_INCLUDE_DIRS)
    target_include_directories(PElib PRIVATE ${PARENT_INCLUDE_DIRS})
endif()

if(DEFINED PARENT_MODULE_DIRS)
    target_include_directories(PElib PRIVATE ${PARENT_MODULE_DIRS})
endif()

install(TARGETS PElib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.mod")

