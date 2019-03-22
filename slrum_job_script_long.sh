#!/bin/bash
#SBATCH -p normal
#SBATCH --nodes=4
# SBATCH --constraint=infiniband,avx2
#SBATCH --tasks-per-node 12
#SBATCH --time=72:00:00

echo "********** CPU-INFO **********"
lscpu
echo "********** XXXXXXXX **********"
# echo "********** Listing TMPDIR **********"

# ls -l
export IMAGE_STORE='/home/alaik/singularityImages'
export UW_ENABLE_TIMING='1'
export IMAGE_NAME='uwgeodynamics-badlands'



echo "********** Run Started **********"

srun -n 48 singularity exec --pwd $PWD $IMAGE_STORE/$IMAGE_NAME.simg  python wedge.py

echo "********** XXXXXXXXXXX **********"

wait
