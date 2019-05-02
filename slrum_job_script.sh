#!/bin/bash
#SBATCH -p short
#SBATCH --nodes=1
# SBATCH --constraint=infiniband,avx2
#SBATCH --tasks-per-node 12
#SBATCH --time=00:05:00

echo "********** CPU-INFO **********"
lscpu
echo "********** XXXXXXXX **********"
# echo "********** Listing TMPDIR **********"

# ls -l
export IMAGE_STORE='/home/alaik/singularityImages'
export UW_ENABLE_TIMING='1'
export IMAGE_NAME='uwgeobadlands-latest'


echo "********** Run Started **********"

srun -n 12 singularity exec --pwd $PWD $IMAGE_STORE/$IMAGE_NAME.simg  python wedge.py

echo "********** XXXXXXXXXXX **********"

wait
