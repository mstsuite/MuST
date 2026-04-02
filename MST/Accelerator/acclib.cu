#include <unistd.h>
#include <iostream>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <nvml.h>
#include "acclib.hpp"
#include "math.h"

void check_cuda_error(cudaError_t result, char const *const func, const char *const file, int const line) {
   if (result != cudaSuccess) {
      fprintf(stderr, "\nCUDA Runtime Error at %s:%d code=%d(%s) \"%s\"\n", file, line,
                      static_cast<unsigned int>(result), cudaGetErrorName(result), func);
      exit(EXIT_FAILURE);
   }
}

void check_cublas_error(cublasStatus_t status, char const *const func, const char *const file, int const line) {
    if (status != CUBLAS_STATUS_SUCCESS) {
       fprintf(stderr, "\ncuBLAS Runtime Error at %s:%d code=%d(%s) \"%s\"\n", file, line,
                       static_cast<unsigned int>(status), cublasGetStatusString(status), func);
       exit(EXIT_FAILURE);
    }
}

extern "C"
void get_node_resources_(int *my_pe, int *num_cpu_cores, int *num_gpu_cards, int *mem_gb){
    *num_cpu_cores = sysconf(_SC_NPROCESSORS_ONLN);    // available cores
    // *num_cpu_cores = sysconf(_SC_NPROCESSORS_CONF); // configured cores

    int deviceCount = 0;
    checkCudaErrors(cudaGetDeviceCount(&deviceCount));
    *num_gpu_cards = deviceCount;

    double mem_min = 1024.0; // Set the minimum of the GPU memory to 1024 GB to begin with
    for (int dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        checkCudaErrors(cudaGetDeviceProperties(&deviceProp, dev)); // Get properties for the current device
        size_t total_bytes = deviceProp.totalGlobalMem;
        double total_mb = (double)total_bytes/(1024.0*1024.0);
        double total_gb = (double)total_bytes/(1024.0*1024.0*1024.0);
        if (*my_pe == 0) {
           std::cout << "\n--- Device " << dev << " ---" << std::endl;
           std::cout << "Name: " << deviceProp.name << std::endl; // The GPU model name
           std::cout << "Compute Capability: " << deviceProp.major << "." << deviceProp.minor << std::endl;
           std::cout << "Total Global Memory: " << total_mb << " MBytes\n" << std::endl;
        }
        mem_min = fmin(mem_min,total_gb);
    }
   if (deviceCount > 0) {
      *mem_gb = round(mem_min+0.005);
   }
   else {
      *mem_gb = 0;
   }
}

nvmlDevice_t device;
unsigned long long start_energy;
double prev_energy;

extern "C"
void initialize_energy_benchmark_() {
    nvmlInit();
    nvmlDeviceGetHandleByIndex(0, &device);

    // Make sure to synchronize to capture the full kernel load
    cudaDeviceSynchronize(); 

    // Get baseline reading before the kernel
    nvmlReturn_t result = nvmlDeviceGetTotalEnergyConsumption(device, &start_energy);
    if (NVML_SUCCESS != result) {
      std::cerr << "Failed to query energy consumption: " << nvmlErrorString(result) << std::endl;
    }
    prev_energy = 0.0e0;
}

extern "C"
void measure_energy_benchmark_(double *present_energy_joules, double *cumulative_energy_joules) {
    // Make sure to synchronize to capture the full kernel load
    cudaDeviceSynchronize(); 

    unsigned long long end_energy;
    // Get reading after completion. Returns end_energy is milliJoules (mJ)
    nvmlReturn_t result = nvmlDeviceGetTotalEnergyConsumption(device, &end_energy);
    if (NVML_SUCCESS != result) {
      std::cerr << "Failed to query energy consumption: " << nvmlErrorString(result) << std::endl;
    }

    // Result in Joules
    double energy_joules = (end_energy - start_energy) / 1000.0;
    *cumulative_energy_joules = energy_joules;
    *present_energy_joules = energy_joules - prev_energy;
    prev_energy = energy_joules;
    
    // std::cout << "Energy consumed by kernels: " << energy_joules << " J" << std::endl;

    /*
    unsigned int power; // Get handle for the first device (index 0)
    // Returns power usage in milliwatts (mW)
    nvmlReturn_t result = nvmlDeviceGetPowerUsage(device, &power);
    if (NVML_SUCCESS == result) {
       std::cout << "Power usage: " << power / 1000.0 << " Watts" << std::endl;
    }
    else {
       std::cerr << "Failed to query power: " << nvmlErrorString(result) << std::endl;
    }
    */
}

extern "C"
void finalize_energy_benchmark_() {
    nvmlShutdown();
}
