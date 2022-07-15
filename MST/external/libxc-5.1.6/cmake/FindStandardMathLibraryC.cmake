# * downloaded Nov 2016 from https://android.googlesource.com/platform/external/eigen/+/master/cmake/FindStandardMathLibrary.cmake
# * changed CXX to C
# * note that full path to libm *not* detected

# SL 2020-10-23: Added more necessary functions for libxc, and changed to a version using standard input

# - Try to find how to link to the standard math library, if anything at all is needed to do.
# On most platforms this is automatic, but for example it's not automatic on QNX.
#
# Once done this will define
#
#  STANDARD_MATH_LIBRARY_FOUND - we found how to successfully link to the standard math library
#  STANDARD_MATH_LIBRARY - the name of the standard library that one has to link to.
#                            -- this will be left empty if it's automatic (most platforms).
#                            -- this will be set to "m" on platforms where one must explicitly
#                               pass the "-lm" linker flag.
#
# Copyright (c) 2010 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.
include(CheckCSourceCompiles)
# a little test program for c++ math functions.
# notice the std:: is required on some platforms such as QNX

# The test program needs to read in the arguments from standard input,
# since otherwise the compiler might optimize away the calls altogether!
set(find_standard_math_library_test_program
"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
int main(int argc, char **argv) {
  double x=atof(argv[1]);
  double y=atof(argv[2]);
  printf(\"sinh(x)=% e\",sinh(x));
  printf(\"atan2(x,y)=% e\",atan2(x,y));
  printf(\"tanh(x)=% e\",tanh(x));
  printf(\"cosh(x)=% e\",cosh(x));
  printf(\"sin(x)=% e\",sin(x));
  printf(\"atan(x)=% e\",atan(x));
  printf(\"exp(x)=% e\",exp(x));
  printf(\"cos(y)=% e\",cos(y));
  printf(\"log(x)=% e\",log(x));
  printf(\"pow(x,y)=% e\",pow(x,y));
  printf(\"round(x)=% e\",round(x));
  printf(\"sqrt(x)=% e\",sqrt(x));
  printf(\"erf(x)=% e\",erf(x));
  return 0;
}")
# C++ test program
# "#include<cmath>
# int main() { std::sin(0.0); std::log(0.0f); }")
# first try compiling/linking the test program without any linker flags
set(CMAKE_REQUIRED_FLAGS "")
set(CMAKE_REQUIRED_LIBRARIES "")
CHECK_C_SOURCE_COMPILES(
  "${find_standard_math_library_test_program}"
  standard_math_library_linked_to_automatically
)
if(standard_math_library_linked_to_automatically)
  message("Test program linked without flags")
  # the test program linked successfully without any linker flag.
  set(STANDARD_MATH_LIBRARY "")
  set(STANDARD_MATH_LIBRARY_FOUND TRUE)
else()
  # the test program did not link successfully without any linker flag.
  # This is a very uncommon case that so far we only saw on QNX. The next try is the
  # standard name 'm' for the standard math library.
  set(CMAKE_REQUIRED_LIBRARIES "m")
  CHECK_C_SOURCE_COMPILES(
    "${find_standard_math_library_test_program}"
    standard_math_library_linked_to_as_m)
  if(standard_math_library_linked_to_as_m)
    # the test program linked successfully when linking to the 'm' library
    set(STANDARD_MATH_LIBRARY "m")
    set(STANDARD_MATH_LIBRARY_FOUND TRUE)

    message("Test program linked with -m")
  else()
    # the test program still doesn't link successfully

    message("Test program did not link with -m")
    set(STANDARD_MATH_LIBRARY_FOUND FALSE)
  endif()
endif()
