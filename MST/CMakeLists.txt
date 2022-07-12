################
# Iman Jehanzeb
################

cmake_minimum_required(VERSION 3.18)

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
message(STATUS "CMake ${CMAKE_VERSION} with max. ${Ncpu} threads")

project(MST LANGUAGES Fortran C)

# make CUDA optional
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
else()
        message(STATUS "No CUDA compiler found")
endif()

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
target_link_libraries(mst_main PUBLIC gfortran)


# Adding Subdirectories
add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(plib)
add_subdirectory(iolib)



# Find MPI Library
find_package(MPI REQUIRED COMPONENTS CXX OPTIONAL_COMPONENTS)
target_link_libraries(mst_main PUBLIC MPI)


# Find OpenMP library
option(BUILD_WITH_OPENMP "Build with OpenMP" ON)
if (BUILD_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(mst_main PUBLIC OpenMP::OpenMP_FORTRAN)
endif ()



# Find MPI library
find_package(MPI OPTIONAL COMPONENTS Fortran OPTIONAL_COMPONENTS)
target_link_libraries(mst_main PUBLIC MPI::MPI_FORTRAN)




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
