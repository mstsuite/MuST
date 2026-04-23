Architecture File Structure
===========================

Architecture files follow a modular structure that defines how the code is built
across different systems (CPU, GPU, compiler toolchains, and libraries).

They are typically organized into the following sections:

---------------------------------------------------------------------

1. Acceleration Control
---------------------------------------

Controls whether GPU acceleration is enabled.

::

    Acceleration = 0   # CPU-only build
    Acceleration = 1   # GPU-enabled build

When enabled, CUDA-related variables and libraries must be defined.

---------------------------------------------------------------------

2. Library Paths
---------------------------------------

Defines paths to external dependencies.

Common variables:
::

    HDF5_PATH
    LIBXC_PATH
    FFTW_PATH
    P3DFFT_PATH
    LUA_PATH
    ACCEL_PATH
    MATH_PATH

Notes:
- Empty paths indicate the package will be built internally.
- External installations should be explicitly provided.
- GPU builds require ACCEL_PATH (CUDA).

---------------------------------------------------------------------

3. Library Linking (LIBS)
---------------------------------------

Specifies libraries required during linking.

Since core external libraries (HDF5, LIBXC, FFTW, P3DFFT) are built internally,
they should NOT be included here.

Typical patterns:

**Minimal CPU (Generic BLAS/LAPACK):**
::

    LIBS += -lblas -llapack -lm -lpthread

**Intel MKL (Sequential):**
::

    LIBS += -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
            -lpthread -lm -ldl

**Intel MKL (MPI + ScaLAPACK):**
::

    LIBS += -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential \
            -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
            -lpthread -lm -ldl

**Intel OneAPI (Simplified):**
::

    LIBS += -qmkl=cluster -lifcore

**GNU + OpenBLAS:**
::

    LIBS += -lopenblas -lm -lpthread

**NVHPC with NVPL:**
::

    LIBS += -Mnvpl -Mscalapack

**NVHPC + explicit BLAS/LAPACK:**
::

    LIBS += -lblas -llapack -lscalapack

**CUDA (Basic GPU support):**
::

    LIBS += -L${NVHPC_ROOT}/cuda/lib64 -lcudart -lcuda \
            -L${NVHPC_ROOT}/math_libs/lib64 -lcublas -lcusolver \
            -lstdc++ -lm

**CUDA + NVHPC (optimized GPU build):**
::

    LIBS += -L${NVHPC_ROOT}/compilers/lib -lblas -llapack \
            -L${NVHPC_ROOT}/cuda/lib64 \
            -lcudart -lnvtx3interop \
            -cuda -cudalib=cublas,cusolver \
            -gpu=cc90,cuda12.9,lineinfo \
            -lstdc++

**MKL + CUDA hybrid:**
::

    LIBS += -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential \
            -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
            -L${NVHPC_ROOT}/cuda/lib64 -lcudart -lcuda \
            -L${NVHPC_ROOT}/math_libs/lib64 -lcublas -lcusolver \
            -lstdc++ -lpthread -lm -ldl

**MKL + CUDA + GNU Fortran runtime:**
::

    LIBS +=  -L${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential \
            -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
            -L${NVHPC_ROOT}/cuda/lib64 -lcudart -lcuda \
            -L${NVHPC_ROOT}/math_libs/lib64 -lcublas -lcusolver \
            -lstdc++ -lpthread -lm -ldl -lgfortran

Notes:
- Choose only ONE math backend (MKL, OpenBLAS, or NVPL).
- GPU builds must include CUDA runtime libraries.
- Add ``-lgfortran`` when mixing GNU Fortran with C/C++ linking.

---------------------------------------------------------------------

4. Compiler Toolchain
---------------------------------------

Defines MPI-enabled compilers.

::

    CC        = mpicc
    CXX       = mpicxx
    F77       = mpif90
    FC        = mpif90
    MPICC     = mpicc
    ACCEL_CXX = nvcc -arch=sm_XX
    ARCHV     = ar -r

Notes:
- Intel toolchains may use mpiicx / mpiifort.
- GPU builds require nvcc.

---------------------------------------------------------------------

5. Compilation and Preprocessor Flags
---------------------------------------

Controls optimization, parallelism, and preprocessing.

::

    FPPDEFS   = -cpp
    CPPDEFS   = -cpp
    FPPFLAGS  = -DMPI -DMaxOutProcs=1

    CFLAGS    = -O3
    CXXFLAGS  = -O3 -std=c++14
    FFLAGS    = -O3 -I.

Optional flags:
::

    -DUSE_SCALAPACK
    -march=native
    -fallow-argument-mismatch

---------------------------------------------------------------------

6. OpenMP Support
---------------------------------------

Defines shared-memory parallelism.

::

    OPT_OPENMP = -fopenmp     # GNU / Intel
    OPT_OPENMP = -mp          # NVHPC

Used in:
::

    CFLAGS, CXXFLAGS, FFLAGS, LD_FLAGS

---------------------------------------------------------------------

7. Linker Settings
---------------------------------------

::

    LD_FLAGS = $(OPT_OPENMP)
    LD       = $(FC) $(LD_FLAGS)

Notes:
- Some builds omit LD_FLAGS.
- GPU builds may include additional CUDA linking flags.

---------------------------------------------------------------------

8. External Package Configuration
---------------------------------------

Used for building dependencies.

::

    HDF5_CONFIG_FLAGS
    LIBXC_CONFIG_FLAGS
    P3DFFT_CONFIG_FLAGS
    FFTW_CONFIG_FLAGS

Typical pattern:
::

    HDF5_CONFIG_FLAGS = --enable-fortran CC=$(CC) FC=$(FC)

Notes:
- MPI and OpenMP are usually enabled.
- Compiler wrappers must match the toolchain.

---------------------------------------------------------------------

9. GPU-Specific Configuration
---------------------------------------

Only required when acceleration is enabled.

::

    ACCEL       = CUDA
    ACCEL_PATH  = /path/to/cuda
    ACCEL_CXX   = nvcc -arch=sm_XX

Typical GPU libraries:
::

    -lcudart -lcublas -lcusolver

Architecture examples:
::

    sm_60, sm_70, sm_80, sm_90
