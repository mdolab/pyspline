cmake_minimum_required(VERSION 3.12)

project(pyspline VERSION 1.5.0)

# Prevent in source builds
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)

if (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
     "Build types are: Debug Release RelWithDebInfo MinSizeRel."
     FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(Fortran)
include(CTest)

# set up directory variables
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(TOPDIR ${CMAKE_SOURCE_DIR})
set(SRC_DIR ${TOPDIR}/src)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Find python and the packages we need
find_package(Python3 3.8 REQUIRED
	COMPONENTS Interpreter Development NumPy)

# Get the variables from python f2py headers
execute_process(
  COMMAND "${Python3_EXECUTABLE}"
  "${SRC_DIR}/f2py/get_f2py.py"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(F2PY_INCLUDE_DIR ${F2PY_INCLUDE_DIR}/src)


# set the include directories
include_directories(
  ${Python3_INCLUDE_DIRS}
  ${Python3_NumPy_INCLUDE_DIRS}
  ${F2PY_INCLUDE_DIR}
  )

message(STATUS "Python 3 Include Dir:    ${Python3_INCLUDE_DIRS}")
message(STATUS "f2py include directory:  ${F2PY_INCLUDE_DIR}")
message(STATUS "Numpy include directory: ${Python3_NumPy_INCLUDE_DIRS}")

set(f2py_module_name "libspline")
set(f2py_module_c "${f2py_module_name}module.c")
set(generated_module_file "${f2py_module_name}${Python3_SOABI}")

# set the fortran precision
if (NOT PRECISION_FLAG)
	get_filename_component(COMPILER_BASE ${CMAKE_Fortran_COMPILER} NAME)
	if (${COMPILER_BASE} STREQUAL "gfortran")
		SET(PRECISION_FLAG -fdefault-real-8 CACHE STRING
			"Specifies the flag to promote real to 8 bytes by default."
			FORCE)
	elseif(${COMPILER_BASE} STREQUAL "ifort")
		SET(PRECISION_FLAG -r8 CACHE STRING
			"Specifies the flag to promote real to 8 bytes by default."
			FORCE)
	else()	
		message(WARNING "Could not find an appropriate flag for promoting 
						default real precision to 8 bytes. Only ifort and 
						gfortran are supported. If you are calling a 
						different compiler and would like to specify this 
						use the -DPRECISION_FLAG=[flag] command line 
						option with cmake.")
	endif()
endif(NOT PRECISION_FLAG)

message(STATUS "Using real promotion flag: ${PRECISION_FLAG}")
add_compile_options(${PRECISION_FLAG})

add_subdirectory(src/)
