# Required version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Project name
PROJECT(Gen1Int C Fortran)

# Changes the default build mode
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers" OFF)

SET(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Header files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

# ConfigParentSettings.cmake configures directories where
# to place mod files and where to find header files
# this is important if there are mod and header file
# dependencies between libraries and where the parent code needs
# to enforce a directory structure to compile correctly
INCLUDE(ConfigParentSettings)

# FFLAGS
#SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -DDEBUG")
#SET(CMAKE_Fortran_FLAGS_RELEASE "-O3")

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qzerosize -qextname -qsuppress=cmpmsg -qflag=w:w")
endif()

# Source codes of Gen1Int library
INCLUDE(Sources)

# Creates library
ADD_LIBRARY(gen1int ${Gen1Int_LIB_SRCS})

# hjaaj Aug 2017: test suite does not compile with
# IBM xlf compiler, because two long statement line
# in refs_boys.h for xlf
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES XL)

# Source codes of Fortran test suite
INCLUDE(TestFortran)

# Creates test suite
ADD_EXECUTABLE(test_gen1int ${Gen1Int_TEST_SRCS})

# Linked libraries for the test suite
TARGET_LINK_LIBRARIES(test_gen1int gen1int)

# Set CDash buildname
SET(BUILDNAME
    "BUILDNAME-not-set" 
    CACHE STRING
    "Name of build on the dashboard")
# Sets CTest own timeout
SET(DART_TESTING_TIMEOUT
    "3600" 
    CACHE STRING
    "Set timeout in seconds for every single test")
INCLUDE(Tests)
INCLUDE(CTest)
ENABLE_TESTING()

endif()

# Installation
INSTALL(TARGETS gen1int ARCHIVE DESTINATION lib)
