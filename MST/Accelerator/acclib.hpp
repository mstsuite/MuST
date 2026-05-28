#ifndef ACCLIB_HPP
#define ACCLIB_HPP

// CUDA error-checking macro for reliable memory allocation and transfers

// CUDA runtime error checking
#define checkCudaErrors(call) {                                      \
    cudaError_t err__ = (call);                                      \
    if (err__ != cudaSuccess) {                                      \
        fprintf(stderr,                                              \
                "CUDA error at %s:%d code=%d (%s)\n",                \
                __FILE__, __LINE__,                                  \
                (int)err__, cudaGetErrorString(err__));              \
        exit(EXIT_FAILURE);                                          \
    }                                                                \
}

// cuBLAS error checking
#define checkCublasErrors(stat) {                                    \
    if ((stat) != CUBLAS_STATUS_SUCCESS) {                           \
        fprintf(stderr,                                              \
                "cuBLAS error at %s:%d code=%d\n",                   \
                __FILE__, __LINE__, (int)(stat));                    \
        exit(EXIT_FAILURE);                                          \
    }                                                                \
}

// cuSOLVER error checking
#define checkCusolverErrors(stat) { \
    if ((stat) != CUSOLVER_STATUS_SUCCESS) { \
        fprintf(stderr, "cuSOLVER error at %s:%d code=%d\n", \
                __FILE__, __LINE__, stat); \
        exit(EXIT_FAILURE); \
    } \
}

#endif
