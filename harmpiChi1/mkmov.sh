#!/bin/bash -l        
#PBS -l walltime=24:00:00
### number of process and nodes to use. each node contains up to 24 processes.
#PBS -l nodes=1:ppn=24
#PBS -m abe 
#PBS -M yoavzack@mail.tau.ac.il
#PBS -N TiltWindMv

module load gcc/gcc-6.2.0
### module load mpi/openmpi-1.10.4
module load mpi/openmpi-2.1.3            
module load intel/intel_com_mkl_mpi2017
module load python/python-2.7-miniconda2.7.15
hostname

cd $PBS_O_WORKDIR

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

python harm_script.py 2> log.mkmov