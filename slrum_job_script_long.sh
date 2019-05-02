#!/bin/bash
#SBATCH -p normal
#SBATCH --nodes=1
# SBATCH --constraint=infiniband,avx2
#SBATCH --tasks-per-node 12
#SBATCH --time=92:00:00

echo "********** CPU-INFO **********"
lscpu
echo "********** XXXXXXXX **********"
# echo "********** Listing TMPDIR **********"

# ls -l
# export IMAGE_STORE='/home/alaik/singularityImages'
# export UW_ENABLE_TIMING='1'
# export IMAGE_NAME='uwgeodynamics-badlands'



echo "********** Run Started **********"

srun -n 12 singularity exec --pwd $PWD docker://arijitlaik/uwgeobadlands:latest  python wedge.py

echo "********** XXXXXXXXXXX **********"

wait
