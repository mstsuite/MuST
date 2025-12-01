# CUDA for NVIDIA
# HIP for AMD

if (USE_ACCELERATOR_CUDA_C)

  target_sources(mst_lib 
              PUBLIC
              acclib.cu
              acclib.hpp
              LSMS_Accel.cu
              )
              
endif()
