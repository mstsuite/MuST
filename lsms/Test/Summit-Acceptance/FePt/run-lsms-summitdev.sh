#!/bin/bash
jobname=LSMS3-FePt
nodes=1
ppn=2
num_threads=8
wt=60
let nmpi=$nodes*$ppn

module load gcc/5.4.0 cuda essl cuda hdf5

# LSMS_ROOT=$HOME/ven101gpfs/LSMS_3-gcc-cuda
LSMS_ROOT=$HOME/MST_Pack/LSMS_3_summitdev

# the run will be for  2 atoms (FePt)
# REPEAT=2
# NUM_ATOMS=$(( 2 * $REPEAT * $REPEAT * $REPEAT ))

# sed "s/REPEAT/$REPEAT/" i_lsms_template > i_lsms_$NUM_ATOMS

# Begin LSF directives
#--------------------------------------
cat > _batch.job <<EOF
#BSUB -J $jobname
#BSUB -P STF006SUMMITDEV
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n ${nmpi}
#BSUB -W $wt
#BSUB -R "span[ptile=$ppn]"
#BSUB -x
#BSUB -env "all,JOB_FEATURE=gpumps"

export OMP_NUM_THREADS=${num_threads}
export OMP_WAIT_POLICY=active
export BIND_THREADS=yes


echo "Running lsms for " $NUM_ATOMS " FePt atoms."

date

export LSMS_ROOT=$HOME/MST_Pack/LSMS_3_summitdev

cd $LSMS_ROOT/Test/Summit-Acceptance/FePt
#cd $MEMBERWORK/stf006
date

rm -f w_fept v_fept
mpirun -mca pml_pami_remote_eager_limit 1048576 -aff off -np ${nmpi} set_device_and_bind.sh $LSMS_ROOT/bin/lsms i_lsms_start

if [ -f w_fept ]; then
  mv w_fept v_fept
  mpirun -mca pml_pami_remote_eager_limit 1048576 -aff off -np ${nmpi} set_device_and_bind.sh $LSMS_ROOT/bin/lsms i_lsms_restart
else
  echo "No 'w_fept' file  after first run of lsms! Check for errors."
fi

EOF
#---------------------------------------
bsub < _batch.job
