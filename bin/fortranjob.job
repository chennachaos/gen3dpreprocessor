#!/bin/bash

# Set the name of the job
# (this gets displayed when you get a list of jobs on the cluster)
#SBATCH --job-name="parallelfort"
#SBATCH --output=fortran-partition.out

# Specify the maximum wall clock time your job can use
# (Your job will be killed if it exceeds this)
#SBATCH --time=5:00:00

# Specify the number of cpu cores your job requires
#SBATCH --ntasks=4

# Specify the amount of memory your job needs per cpu-core (in Mb)
# (Your job will be killed if it exceeds this for a significant length of time)
#SBATCH --mem-per-cpu=5000


# Set up the environment
module purge
module load hpcw
module load petsc/3.7.5
#module load compiler/intel/2018/3
#module load mpi/intel/2018/3

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/libraries/impi/5.0.1.035/lib64
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/s.engkadac/mylibs/parmetis-4.0.3-install/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/libraries/petsc/3.7.5/el6/AVX/intel-16.0/intel-5.1/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/libraries/petsc/3.7.5/el6/SSE4.2/intel-16.0/intel-5.1/lib

# Run the application
echo My job is started

mpirun ./partfort m6 10

echo My job has finished





