#!/bin/bash
jobname=LSMS3-Fe
nodes=4
ppn=4
num_threads=4
wt=60
let nmpi=$nodes*$ppn

module load gcc/5.4.0 cuda essl cuda hdf5

# LSMS_ROOT=$HOME/ven101gpfs/LSMS_3-gcc-cuda
LSMS_ROOT=$HOME/MST_Pack/LSMS_3_summitdev

# the run will be for  2 x REPEAT^3  Fe atoms
REPEAT=4
NUM_ATOMS=$(( 2 * $REPEAT * $REPEAT * $REPEAT ))

sed "s/REPEAT/$REPEAT/" i_lsms_template > i_lsms_$NUM_ATOMS

# Begin LSF directives
#--------------------------------------
cat > _batch.job <<EOF
#BSUB -J $jobname
#BSUB -P STF006SUMMITDEV
#BSUB -o %J.out
#BSUB -e %J.err
##BSUB -n ${nmpi}
#BSUB -nnodes ${nodes}
#BSUB -W $wt
##BSUB -R "span[ptile=$ppn]"
##BSUB -x
##BSUB -env "all,JOB_FEATURE=gpumps"

export OMP_NUM_THREADS=${num_threads}
export OMP_WAIT_POLICY=active
export BIND_THREADS=yes

echo "Running lsms for " $NUM_ATOMS " Fe atoms."

date

# mpirun -aff off -np ${nmpi} set_device_and_bind.sh $LSMS_ROOT/bin/lsms i_lsms_$NUM_ATOMS
# mpirun -mca pml_pami_remote_eager_limit 1048576 -aff off -np ${nmpi} set_device_and_bind.sh $LSMS_ROOT/bin/lsms i_lsms_$NUM_ATOMS

jsrun -n${nmpi} -a1 -g1 -r${ppn} -c${num_threads} $LSMS_ROOT/bin/lsms i_lsms_$NUM_ATOMS

EOF
#---------------------------------------
bsub  _batch.job
