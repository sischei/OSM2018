#!/bin/sh
#BATCH --job-name=job1
#SBATCH --output=job1.out
#SBATCH --error=job1.err

#SBATCH --time=00:00:10

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16



module unload openmpi 
module load mpi4py/1.3+python-2.7-2015q2

mpirun python bcast.py

