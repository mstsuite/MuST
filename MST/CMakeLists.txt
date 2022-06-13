################
# Iman Jehanzeb
################

cmake_minimum_required(VERSION 3.18)

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
message(STATUS "CMake ${CMAKE_VERSION} with max. ${Ncpu} threads")

project(MST LANGUAGES Fortran C CXX CUDA)

# Defining mst_main
add_executable(mst_main)

set_target_properties(mst_main PROPERTIES OUTPUT_NAME "MST")

set_target_properties(mst_main PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )


# Defining mst_lib 
add_library(mst_lib)

set_target_properties(mst_lib PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
        
target_include_directories(mst_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(mst_main PUBLIC mst_lib)
target_link_libraries(lsmscore PUBLIC gfortran)


# Adding Subdirectories
add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(plib)
add_subdirectory(iolib)



# Find MPI Library
find_package(MPI REQUIRED COMPONENTS CXX OPTIONAL_COMPONENTS)
target_link_libraries(mst_main PUBLIC MPI::MPI_CXX)


# Find OpenMP library
option(BUILD_WITH_OPENMP "Build with OpenMP" ON)
if (BUILD_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(mst_main PUBLIC OpenMP::OpenMP_CXX)
endif ()


# Needed for gcc 10
target_compile_options(mst_lib PUBLIC "-fallow-argument-mismatch")
        

# Build LibXC 
option(USE_LIBXC "Use the libxc libary" ON)
if (USE_LIBXC)
    message(STATUS "LIBXC is used")
    include(libxc)
else ()
    message(STATUS "LIBXC is not used")
endif ()


# Build mjson 
include(mjson)


# CUDA
if (USE_ACCELERATOR_CUDA_C)

    target_compile_definitions(mst_main PUBLIC ACCELERATOR_CUDA_C)

    option(CUDA_DEBUG "Debug CUDA routines" OFF)

    if (CUDA_DEBUG)
        target_compile_definitions(mst_main PUBLIC "CUDA_DEBUG")
    endif ()

        
# Build with OpenMP
if (BUILD_WITH_OPENMP)
        message(STATUS "OpenMP is also defined for GPU")
        target_compile_options(mst_main PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -fopenmp>)
    endif ()

    if (NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 14)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif ()

    if (NOT DEFINED CMAKE_CUDA_SEPARABLE_COMPILATION)
        set(CMAKE_CUDA_SEPARABLE_COMPILATION OFF)
    endif ()

    find_package(CUDAToolkit REQUIRED)

    message(STATUS "CUDA Version: " ${CUDAToolkit_VERSION})
    message(STATUS "CUDA Include dirs: " ${CUDAToolkit_INCLUDE_DIRS})
    message(STATUS "CUDA Library: " ${CUDAToolkit_LIBRARY_DIR})

    target_link_libraries(mst_main PUBLIC CUDA::toolkit)

    # Check if the necessary libraries have been found
    if (TARGET CUDA::cublas)
        target_link_libraries(mst_main PUBLIC CUDA::cublas)
    else ()
        message(FATAL_ERROR "CuBLAS was not defined or found")
    endif ()

    if (TARGET CUDA::cusolver)
        target_link_libraries(mst_main PUBLIC CUDA::cusolver)
    else ()
        message(FATAL_ERROR "Cusolver was not defined or found")
    endif ()
endif ()


if (DEFINED MST_LINEAR_SOLVER_DEFAULT)
    message(STATUS "Default linear solver: " ${MST_LINEAR_SOLVER_DEFAULT})
    set(MST_LINEAR_SOLVER_DEFAULT ${MST_LINEAR_SOLVER_DEFAULT} CACHE STRING "Default linear solver")
    target_compile_definitions(mst_main PUBLIC "MST_LINEAR_SOLVER_DEFAULT=${MST_LINEAR_SOLVER_DEFAULT}")
endif ()

if (DEFINED MST_BUILD_KKR_MATRIX_DEFAULT)
    message(STATUS "Default kkr matrix build: " ${MST_BUILD_KKR_MATRIX_DEFAULT})
    set(MST_BUILD_KKR_MATRIX_DEFAULT ${MST_BUILD_KKR_MATRIX_DEFAULT} CACHE STRING "Default linear solver")
    target_compile_definitions(mst_main PUBLIC "MST_BUILD_KKR_MATRIX_DEFAULT=${MST_BUILD_KKR_MATRIX_DEFAULT}")
endif ()


option(ARCH_IBM "Is IBM architecture" OFF)

if (ARCH_IBM)
    message(STATUS "Use the IBM architecture")
    target_compile_definitions(mst_main PUBLIC "ARCH_IBM")
endif ()

# Find MPI library
find_package(MPI REQUIRED COMPONENTS CXX OPTIONAL_COMPONENTS)
target_link_libraries(mst_main PUBLIC MPI::MPI_CXX)


# Find HDF5
message(STATUS "Find HDF5")
set(HDF5_PREFER_PARALLEL false)
find_package(HDF5 REQUIRED)

message(STATUS "HDF5 C Libaries: " "${HDF5_C_LIBRARIES}")
message(STATUS "HDF5 LIBRARIES: " "${HDF5_LIBRARIES}")
message(STATUS "HDF5 INCLUDE DIRS: " "${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 Version: " "${HDF5_VERSION}")

target_link_libraries(mst_main PUBLIC HDF5::HDF5)


# Linear algebra libraries
option(USE_ESSL "Use IBM ESSL library" OFF)

if (USE_ESSL)
    find_package(ESSL REQUIRED)
    set(ESSL_WORKAROUND True)

    target_link_libraries(mst_main PUBLIC ESSL)
else ()
    # Lapack
    option(SEARCH_LAPACK "Search for LAPACK" ON)
    if (NOT SEARCH_LAPACK)
        add_library(LAPACK::LAPACK STATIC IMPORTED GLOBAL)
        set_target_properties(LAPACK::LAPACK PROPERTIES IMPORTED_LOCATION ${LAPACK_LIBRARIES})
        if (DEFINED LAPACK_LINKER_FLAGS)
            target_compile_options(LAPACK::LAPACK INTERFACE ${LAPACK_LINKER_FLAGS})
        endif ()
        if (DEFINED LAPACK_INCLUDE_DIR)
            target_include_directories(LAPACK::LAPACK INTERFACE ${LAPACK_INCLUDE_DIR})
        endif ()
    else ()
        if (NOT DEFINED BLA_VENDOR)
            set(BLA_VENDOR "Generic")
        endif ()
        find_package(LAPACK REQUIRED)
    endif ()
    target_link_libraries(mst_main PUBLIC LAPACK::LAPACK)


# BLAS
option(SEARCH_BLAS "Search for BLAS" ON)
if (NOT SEARCH_BLAS)
        add_library(BLAS::BLAS STATIC IMPORTED GLOBAL)
        set_target_properties(BLAS::BLAS PROPERTIES IMPORTED_LOCATION ${BLAS_LIBRARIES})
        if (DEFINED BLAS_LINKER_FLAGS)
                target_compile_options(BLAS::BLAS INTERFACE ${BLAS_LINKER_FLAGS})
        endif ()
        if (DEFINED BLAS_INCLUDE_DIR)
                target_include_directories(BLAS::BLAS INTERFACE ${BLAS_INCLUDE_DIR})
        endif ()

        else ()
                if (NOT DEFINED BLA_VENDOR)
                        set(BLA_VENDOR "Generic")
                endif ()
                find_package(BLAS REQUIRED)
        endif ()
        target_link_libraries(mst_main PUBLIC BLAS::BLAS)
endif ()
