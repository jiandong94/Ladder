#!/bin/bash
#PBS -N ob_UORUO
#PBS -l nodes=1:ppn=8
#PBS -l walltime=200:00:00
#PBS -l mem=20G
#PBS -j oe
#PBS -k oe
#PBS -V

APP=observer
INPUTFILE=IF_U0RU0


# Number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes
#export OMP_NUM_THREADS=8

WORKDIR=$PBS_O_WORKDIR
echo Working directory is $WORKDIR
cd $WORKDIR

date
./$APP $INPUTFILE  2>&1 | tee o${ID}.${PBS_JOBNAME}
date

