# CUDA for NVIDIA
# HIP for AMD

if (USE_ACCELERATOR_CUDA_C)

  target_sources(mst_lib 
              PUBLIC
              inverse.cu
              inverse.h
              DeviceStorage.hpp
              DeviceStorage.cu
              operations.h
              zblock_lu_cuda_c.cu
              cudaDoubleComplex.hpp
              initialize_matinv_cuda.F90
              finalize_matinv_cuda.F90
              cudaCheckError.hpp
              )
              
endif()


#Array3d.hpp
#BLAS.hpp
#Complex.hpp
#DeviceArray3d.hpp
#DeviceMatrix.hpp
#DeviceVector.hpp
#LAPACK.hpp
#Matrix.hpp
#Real.hpp
#TypeTraits.hpp
#finalize_matinv_accel.F90
#initialize_matinv_accel.F90
#interface_dstore.cpp
#lapack.h
#zblock_lu_CPU.F
#zblock_lu_accel.F90
