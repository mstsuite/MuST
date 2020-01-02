#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

void *deviceStorage;
void *allocateDStore(void);
void freeDStore(void * d_store);
int initDStore(void * d_store,int kkrsz_max, int nspin, int numLIZ, int nthreads);

#ifdef MAX_GPU_THREADS
#define _max_gpu_threads MAX_GPU_THREADS
#else
#define _max_gpu_threads 1
#endif

extern "C"
void initialize_dstore_(int nthreads, int *num_blocks, int *max_block_sz) {
    // printf("LIZ = %d\n",num_blocks[0]);
    // printf("KKRSZ*nspin_cant = %d\n",max_block_sz[0]);
    int n = std::min(nthreads,_max_gpu_threads);
    if (n < 1) {
       printf("The number of GPU threads can not be negative...: %d\n",n);
       exit(1);
    }
    else if (initDStore(deviceStorage,max_block_sz[0],1,num_blocks[0],n) != 0) {
       printf("Return error by initDStore...\n");
       exit(1);
    }
}

extern "C"
void finalize_dstore_() {
    freeDStore(deviceStorage);
}

extern "C"
void allocate_dstore_() {
    allocateDStore();
}
