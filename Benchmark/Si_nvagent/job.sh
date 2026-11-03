#!/bin/bash
#SBATCH -J job # Job name
#SBATCH -o job.o%j    # Name of stdout output file
#SBATCH -e job.e%j    # Name of stderr error file
#SBATCH -p gh         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes
#SBATCH -n 64            # Total # of mpi tasks
#SBATCH -t 00:10:00       # Run time (hh:mm:ss)
#SBATCH -A DMR22013          # Project allocation

ml purge
ml use /home1/apps/nvidia/modulefiles
ml nvhpc-hpcx-cuda13/26.1

# Set MPS environment variables

export PATH="/work/08910/nnepal/vista/MuST/bin/try1:$PATH"
#export LD_PRELOAD=/work/08910/nnepal/vista/scilib-accel/scilib-dbi.so:/work/08910/nnepal/vista/GEMMul8/lib/libgemmul8.so

export OMPI_MCA_coll="^hcoll"
export NO_STOP_MESSAGE=1
export UCX_WARN_UNUSED_ENV_VARS=n
export UCX_LOG_LEVEL=ERROR
ibrun mst2 < i_mst
