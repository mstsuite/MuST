
OBJ       = accelerator_initialize.o accelerator_finalize.o accelerator_interface.o Accelerator.o \
            buildKKRMatrix_kernel_gpu.o buildKKRMatrix_gpu.o  zblock_lu_cuda_c.o fortran.o inverse.o DeviceStorage.o 

all: libAccelerator.a

%.o: %.F
	$(F77) -c $(FFLAGS) $< -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $(INC_PATH) $< -o $@

%.o: %.cpp
	$(CXX) $(INC_PATH) -c $(CFLAGS) $< -o $@

buildKKRMatrix_gpu.o: buildKKRMatrix_gpu.cu
ifdef CUDA_CXX
	$(CUDA_CXX) -c buildKKRMatrix_gpu.cu $(INC_PATH) -I. $(OPT)
else
	$(CXX) -c no_gpu_buildKKRMatrix.cpp -o buildKKRMatrix_gpu.o
endif

buildKKRMatrix_kernel_gpu.o: buildKKRMatrix_kernels.cu makebgij_device.hpp plglmax_device.hpp  DeviceMatrix.hpp DeviceVector.hpp DeviceArray3d.hpp no_gpu_buildKKRMatrix.cpp
ifdef CUDA_CXX
	$(CUDA_CXX) -c buildKKRMatrix_kernels.cu $(INC_PATH) -I. $(OPT) -o buildKKRMatrix_kernel_gpu.o
else
	$(CXX) -c no_gpu_kernel_buildKKRMatrix.cpp -o buildKKRMatrix_kernel_gpu.o
endif

zblock_lu_cuda_c.o: zblock_lu_cuda_c.cu 
ifdef CUDA_CXX
	$(CUDA_CXX) -c zblock_lu_cuda_c.cu $(INC_PATH) -I. -I $(MPICH_DIR)/include $(OPT)
else
	touch zblock_lu_cuda_c.o
endif

DeviceStorage.o:  DeviceStorage.cu
ifdef CUDA_CXX
	$(CUDA_CXX) -c DeviceStorage.cu $(INC_PATH) -I. $(OPT)
else
	touch DeviceStorage.o
endif

inverse.o:  inverse.cu inverse.h operations.h
ifdef CUDA_CXX
	$(CUDA_CXX) -c inverse.cu $(INC_PATH) -I. $(OPT)
else
	touch inverse.o
endif


clean:
	rm -f *.o *.a

libAccelerator.a: $(OBJ)
	$(AR) -r libAccelerator.a $(OBJ)
	cp libAccelerator.a $(TOP_DIR)/lib

