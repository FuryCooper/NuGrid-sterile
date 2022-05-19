#!/bin/bash
#PBS -l nodes=16:ppn=7,walltime=12:00:00
#PBS -N 2LPTic
ARCH=$(uname -m)
export P4_GLOBMEMSIZE=40000000

/usr/mpi/gcc/openmpi-1.2.2-1/bin/mpirun -mca btl ib,tcp,self  -hostfile $PBS_NODEFILE -np 112 /home/rs123/2LPT/2lpt.sh /home/rs123/2LPT/2lpt_Carmen.param

exit 0;
