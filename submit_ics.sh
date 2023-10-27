#!/bin/bash
###### Job name ######
#PBS -N music_level11_100
###### Output files ######
#PBS -o s_ics_100.out
#PBS -e s_ics.err
###### Number of nodes and cores ######
#PBS -l nodes=1:ncpus=64
###### Specify how many hours do you need ######
#PBS -l walltime=2:00:00
###### Queue name ######
#PBS -q ct64
###### Sends mail to yourself when the job begins and ends ######
#PBS -M myho@asiaa.sinica.edu.tw
#PBS -m be
###### Sandbox ######
#PBS -W sandbox=PRIVATE
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
# module load gcc/11.2.0 hdf5/1.12.1 intel/2022.1.0 mpich/4.0.1 python/3.10.2 fftw/3.3.10 gsl/2.7
module load gcc/11.2.0 hdf5/1.12.1 intel/2022.0.2 mpi/2021.5.1 fftw/2.1.5 gsl/2.7
# module load gcc/11.2.0 openmpi/4.1.4 szip hdf5/1.12.1 python/3.10.2 fftw/2.1.5 gsl/2.7

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/opt/fftw/2.1.5/mpi/2021.5.1/intel/2022.0.2/lib

###### Run your jobs with parameters ######
if [ -n "$PBS_NODEFILE" ]; then
  if [ -f $PBS_NODEFILE ]; then
    NPROCS=`wc -l < $PBS_NODEFILE`
  fi
fi

#$OPENMPI_HOME/bin/mpirun -v -machinefile $PBS_NODEFILE -np $NPROCS ./a.out foo=bar

HDF5_DISABLE_VERSION_CHECK=1 mpirun -np 64 ./GIZMO music/params_music_100.txt
