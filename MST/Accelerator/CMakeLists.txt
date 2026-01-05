# CUDA for NVIDIA
# HIP for AMD

if (USE_ACCELERATOR_CUDA_C)

  target_sources(mst_lib 
              PUBLIC
              acclib.cu
              acclib.hpp
              accmath.hpp  
              util_mat.cu
              LSMS_Accel.cu
              cusolver_KKR_c.cu  
              invertMatrixKKR_CUDA.F90 
              )
              
endif()
