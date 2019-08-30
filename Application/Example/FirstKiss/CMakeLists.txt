### CMake project file for Example
cmake_minimum_required (VERSION 2.8.7)


### set default compiler as intel compiler for MPI if exist
if(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_Fortran_COMPILER NAMES mpiifort DOC "Fortran compiler")
  if(CMAKE_Fortran_COMPILER)
    message(STATUS "Intel Fortran Compiler ${CMAKE_Fortran_COMPILER} is found")
  endif(CMAKE_Fortran_COMPILER)
endif(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
if(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_C_COMPILER NAMES mpiicc DOC "C compiler")
  if(CMAKE_C_COMPILER)
    message(STATUS "Intel C Compiler ${CMAKE_C_COMPILER} is found")
  endif(CMAKE_C_COMPILER)
endif(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
if(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_CXX_COMPILER NAMES mpiicpc DOC "CXX compiler")
  if(CMAKE_CXX_COMPILER)
    message(STATUS "Intel CXX Compiler ${CMAKE_CXX_COMPILER} is found")
  endif(CMAKE_CXX_COMPILER)
endif(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)

### set default compiler as intel compiler if exist
if(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_Fortran_COMPILER NAMES ifort DOC "Fortran compiler")
  if(CMAKE_Fortran_COMPILER)
    message(STATUS "Intel Fortran Compiler ${CMAKE_Fortran_COMPILER} is found")
  endif(CMAKE_Fortran_COMPILER)
endif(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
if(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_C_COMPILER NAMES icc DOC "C compiler")
  if(CMAKE_C_COMPILER)
    message(STATUS "Intel C Compiler ${CMAKE_C_COMPILER} is found")
  endif(CMAKE_C_COMPILER)
endif(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
if(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_CXX_COMPILER NAMES icpc DOC "CXX compiler")
  if(CMAKE_CXX_COMPILER)
    message(STATUS "Intel CXX Compiler ${CMAKE_CXX_COMPILER} is found")
  endif(CMAKE_CXX_COMPILER)
endif(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)


### set project name
project (Example)
enable_language (Fortran)
enable_language (C)


### set SHOWERMCTOP
set (SHOWERMCTOP $ENV{SHOWERMCTOP})
if (NOT SHOWERMCTOP)
  message(FATAL_ERROR "$SHOWERMCTOP is not defined")
endif (NOT SHOWERMCTOP)

### set COSMOSTOP AND COSMOSINC
set (COSMOSTOP $ENV{COSMOSTOP})
if (NOT COSMOSTOP)
  message(FATAL_ERROR "$COSMOSTOP is not defined")
endif (NOT COSMOSTOP)
set (COSMOSINC ${COSMOSTOP}/cosmos)

### set LIBLOFT AND LIBLOFTINC
set (LIBLOFT $ENV{LIBLOFT})
if (NOT LIBLOFT)
  message(FATAL_ERROR "$LIBLOFT is not defined")
endif (NOT LIBLOFT)
set (LIBLOFTINC ${LIBLOFT}/Header)


### make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)


### default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
  "Choose the installation directory; by default it installs in the COSMOS directory."
  FORCE)


### FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


### set compile option
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (ARCH LinuxGfort)
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -std=legacy -fno-align-commons -O3 -D${ARCH}")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -std=legacy -fno-align-commons -O0 -D${ARCH} -g")
  set (CMAKE_C_FLAGS_RELEASE "-O3 -D${ARCH}")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -D${ARCH} -g")
elseif (Fortran_COMPILER_NAME MATCHES "mpiifort.*")
  # mpiifort
  set (ARCH MPIIFC64)
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O3 -D${ARCH} -DDOMPI")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O0 -D${ARCH} -g ")
  set (CMAKE_C_FLAGS_RELEASE "-O3 -D${ARCH}")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -D${ARCH} -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (ARCH PCLinuxIFC64)
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O3 -D${ARCH}")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O0 -D${ARCH} -g ")
  set (CMAKE_C_FLAGS_RELEASE "-O3 -D${ARCH}")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -D${ARCH} -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


### dump variables
message (STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message (STATUS "C compiler: " ${CMAKE_C_COMPILER})
message (STATUS "C++ compiler: " ${CMAKE_CXX_COMPILER})
message (STATUS "$SHOWERMCTOP is ${SHOWERMCTOP}")
message (STATUS "$LIBLOFT is ${LIBLOFT}")
message (STATUS "$COSMOSTOP is ${COSMOSTOP}")
message (STATUS "$ARCH is ${ARCH}")


include_directories(${LIBLOFTINC})
include_directories(${COSMOSINC})
set (CMAKE_Fortran_MODULE_DIRECTORY "${SHOWERMCTOP}/lib/${ARCH}")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${SHOWERMCTOP}/lib/${ARCH}")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${SHOWERMCTOP}/lib/${ARCH}")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set (LIBSHOWERMC showermc)
link_directories("${SHOWERMCTOP}/lib/${ARCH}")

add_executable(cosmos${ARCH} chook.f)
target_link_libraries(cosmos${ARCH} ${LIBSHOWERMC})
